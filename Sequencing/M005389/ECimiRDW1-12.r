require("pacman")

# Load the packages
p_load(
  data.table, scales, edgeR, statmod, poolr,
  pheatmap, svglite, ggplot2, ggrepel, RColorBrewer, tidyverse, magrittr, ggpubr, ggallin
)

show_col(brewer.pal(12, "Paired"))

# Load the data, first column is sample, second column is spacer, third column is count
data <- fread(
  "Sequencing/M005389/ECimiRDW1-12.tsv.gz",
  header = FALSE, sep = "\t", col.names = c("sample", "spacer", "count")
) %>% filter(!(spacer %like% "\\*"))

freezer_stock <- fread(
  "Sequencing/Freezer_Stocks/ECO_freezer.tsv.gz",
  header = TRUE, sep = "\t", col.names = c("sample", "spacer", "count")
) %>% filter(!(spacer %like% "\\*"))


# load the design
design <- fread(
  "Sequencing/M005389/design.tsv",
  header = TRUE, sep = "\t"
)

# inner_join the design with the data
data <- data %>% inner_join(design, by = "sample")

# load the targets
targets <- fread(
  "Organisms/A_E_coli.targets.tsv.gz",
  header = TRUE, sep = "\t"
) %>%
  mutate(
    type = case_when(
      mismatches == 0 ~ "perfect",
      mismatches == 1 ~ "mismatch",
      mismatches == "None" ~ "control"
    )
  )

targets[type != "control", target := toupper(target)]

# targets[sp_dir == tar_dir, type := "control"]
# targets[overlap < 20, type := "control"]

essentials <- targets %>% filter(type == "mismatch") %>% pull(locus_tag) %>% unique()

targets[locus_tag %in% essentials & type == "perfect", type := "perfect essential"]
targets[spacer %in% targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]

# create a data.table with the types for each spacer, i.e., perfect, mismatch, control
types <- targets %>%
  select(spacer, chr, target, mismatches, type) %>%
  unique()

data <- data %>% inner_join(types, by = c("spacer"))

# draw a log CPM replicate correlation plot
replicate_plot <- data %>%
  group_by(sample) %>%
  mutate(cpm = 1e6 * count / sum(count)) %>%
  data.table() %>%
  dcast(type + spacer + imipenem + induced ~ replicate, value.var = "cpm", fill = 0) %>%
  arrange(desc(type)) %>%
  ggplot(aes(x = A, y = B)) +
  geom_point(aes(color = type), size = 1, alpha = 0.5) +
  scale_color_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  stat_cor(method = "pearson", aes(label = after_stat(label)), label.x = 0, label.y = 5, na.rm = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(imipenem ~ induced, labeller = labeller(
    .rows = function(rows) paste("Imipenem:", rows),
    .cols = function(cols) ifelse(cols == "TRUE", "Induced", "Not induced")
  )) +
  theme_bw() +
  labs(color = "Type") +
  theme(strip.text.y = element_text(angle = 0))

print(replicate_plot)

# 0.25 looks bad, so we'll remove it
data <- data %>% filter(!(imipenem == 0.25))

# fill in missing values with 0
full_data <- data %>%
  select(sample, spacer, count) %>%
  dcast(spacer ~ sample, value.var = "count", fill = 0) %>%
  melt(value.name = "count", variable.name = "sample", id.vars = "spacer")

# inner_join with design and targets
full_data <- full_data %>%
  inner_join(design, by = "sample") %>%
  inner_join(types, by = "spacer")

# create names for the samples based on the design
design_names <- full_data %>%
  select(sample, induced, imipenem, replicate, doublings) %>%
  unique() %>%
  mutate(induction = ifelse(induced == "TRUE", "induced", "uninduced")) %>%
  arrange(imipenem, desc(induction), replicate) %>%
  mutate(
    verbose = paste(induction, imipenem, replicate),
    group = paste(induction, imipenem, sep = "_")
  ) %>%
  mutate(group = gsub("\\.", "_", group)) %>%
  select(-induction)

# create factor for the design
design_names$group <- factor(design_names$group, levels = unique(design_names$group))

# explicitly reorder design_names by the groups (probably not necessary, but just in case)
design_names <- design_names %>% arrange(group)

# then factor the design_names$sample by existing order, i.e., the order of the design_names$group
design_names$sample <- factor(design_names$sample, levels = design_names$sample)

# create design matrix for edgeR with the groups
design_matrix <- model.matrix(~ 0 + design_names$group) %>%
  set_colnames(levels(design_names$group)) %>%
  set_rownames(design_names$sample)

# explicitly reorder the full_data by the design_groups
full_data <- design_names %>%
  inner_join(full_data) %>%
  arrange(group)

# list of spacer names
spacers <- full_data %>%
  select(spacer) %>%
  unique() %>%
  pull()

# create a matrix of counts for edgeR
full_data_matrix <- full_data %>%
  dcast(spacer ~ factor(sample, levels = unique(sample)), value.var = "count") %>%
  select(-spacer) %>%
  data.matrix() %>%
  set_rownames(spacers)

# create a DGEList object
dge <- DGEList(counts = full_data_matrix, group = design_names$group)

# filterByExpr
keep <- dge %>%
  filterByExpr(design = design_matrix, group = design_names$group)

# keep only the spacers that pass the filterByExpr
dge <- dge[keep, , keep.lib.sizes = TRUE]

# # normalize
dge <- calcNormFactors(dge)

# scaleOffset using the doublings
# dge <- scaleOffset(dge, design_names$doublings)


# estimate dispersion
dge <- estimateDisp(dge, design_matrix)

# create a contrast matrix
contrasts <- makeContrasts(
  induction_only = induced_0 - uninduced_0,
  uninduced_drift = uninduced_0_125 - uninduced_0,
  imipenem_full = induced_0_125 - uninduced_0_125,
  imipenem_partial = induced_0_125 - induced_0,
  levels = design_matrix
)

# fit the glmQLFTest
fit <- glmQLFit(dge, design = design_matrix, robust = TRUE)

# create a single data.table with all the results, going through each contrast, one at a time
results <- lapply(colnames(contrasts), function(contrast) {
  topTags(glmQLFTest(fit, contrast = contrasts[, contrast]), n = Inf) %>%
    use_series(table) %>%
    data.table(keep.rownames = "spacer") %>%
    inner_join(targets) %>%
    mutate(contrast = contrast)
}) %>%
  rbindlist()

# add the contrast as a factor, so that the order is preserved
results$contrast <- factor(results$contrast, levels = colnames(contrasts))

# draw volcano plots with facets
volcano_plots <- results %>%
  arrange(desc(type)) %>%
  ggplot(aes(x = logFC, y = FDR)) +
  scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
  geom_point(aes(color = type), size = 2, alpha = 0.5) +
  scale_color_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C")) +
  theme_bw() +
  labs(color = "Type") +
  facet_grid(contrast ~ ., scales = "free") +
  theme(strip.text.y = element_text(angle = 0))

print(volcano_plots)

# create gene-level summary for perfect spacers for each contrast

perfect_median_results <- results %>%
  filter(type %in% c("perfect", "perfect essential") & locus_tag != "None") %>%
  group_by(contrast, locus_tag, type) %>%
  summarize(
    logFC = median(logFC),
    FDR = poolr::stouffer(FDR)$p,
    .groups = "drop"
  )

# create faceted volcano plots for perfect median results
volcano_plots <- perfect_median_results %>%
  ggplot(aes(x = logFC, y = FDR)) +
  scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
  geom_point(size = 2, alpha = 0.5, aes(color = type)) +
    scale_color_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C")) +
  theme_bw() +
  facet_grid(contrast ~ ., scales = "free") +
  theme(strip.text.y = element_text(angle = 0))

print(volcano_plots)

definitions <- fread("Organisms/E_coli_genes_from_string.tsv",
  header = TRUE,
  sep = "\t",
  col.names = c("string_id", "gene_name", "size", "annotation")
) %>%
  mutate(locus_tag = gsub(".*\\.", "", string_id))

results_summary <- results %>%
  filter(type %in% c("perfect", "perfect essential")) %>%
  dcast(locus_tag ~ contrast, value.var = "logFC", fun.aggregate = median) %>%
  left_join(definitions)


# results %>%
#       filter(type == "perfect essential") %>%
#       select(target, logFC) %>%
#       rename(perfect_LFC = logFC) %>%
#       inner_join(results %>% filter(type == "mismatch")) %>% rename(mismatch_LFC = logFC) %>%
#       ggplot(aes(x = perfect_LFC, y = mismatch_LFC)) +
#       geom_point(aes(colour = as.numeric(pmin(y_pred, 1))), alpha = 0.1) +
#       scale_color_gradientn(colors = viridis::magma(100)) +
#       facet_wrap(~contrast)