require("pacman")

# Load the packages
p_load(
  data.table, scales, edgeR, statmod, poolr, ggtext,
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
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_point(aes(color = type), size = 1, alpha = 0.5, shape = 16) +
  scale_color_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C")) +
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


# full_data <- full_data %>% inner_join(median_spacers)

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
dge <- estimateGLMRobustDisp(dge, design_matrix)

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


definitions <- fread("Organisms/E_coli_genes_from_string.tsv",
  header = TRUE,
  sep = "\t",
  col.names = c("string_id", "gene_name", "size", "annotation")
) %>%
  mutate(locus_tag = gsub(".*\\.", "", string_id))

results_summary <- results %>%
  filter(target != "None") %>%
  # filter(type %in% c("perfect", "perfect essential")) %>%
  dcast(locus_tag + type ~ contrast, value.var = "logFC", fun.aggregate = median) %>%
  left_join(definitions %>% select(-annotation))

median_spacers <- data %>% select(sample, spacer, count) %>%
  # rbind(freezer_stock) %>% inner_join(targets %>% select(spacer, target) %>% unique()) %>% 
  inner_join(types) %>% 
  # filter(type %in% c("mismatch", "perfect")) %>%
  group_by(sample) %>%
  mutate(cpm = cpm(count), lcpm = log(cpm)) %>% 
  group_by(target) %>%
  mutate(LMT_count = log(mean(cpm)), t_deviance = lcpm - LMT_count, t_sd = sd(cpm)) %>% 
  filter(abs(t_deviance) == min(abs(t_deviance)) & t_sd == min(t_sd)) %>%
  inner_join(targets %>% filter(sp_dir != tar_dir & overlap == 20 & offset >= 0)) %>%
  group_by(locus_tag) %>%
  mutate(LMG_count = log(mean(cpm)), g_deviance = t_deviance - LMG_count, g_sd = sd(cpm)) %>%
  filter(abs(g_deviance) == min(abs(g_deviance)) & g_sd == min(g_sd)) %>%
  ungroup %>%
  select(locus_tag, spacer) %>% unique()



overall_median_results <- results %>% 
      inner_join(median_spacers) %>%
      filter(target != "None") %>% data.table() %>%
      dcast(locus_tag + type ~ contrast, value.var = "logFC", fun.aggregate = median) %>%
      left_join(definitions %>% select(-annotation)) %>% 
      arrange(imipenem_partial)

# create faceted volcano plots for perfect median results
volcano_plots <- median_spacers %>%
  inner_join(results) %>%
  mutate(FDR = ifelse(FDR == 1, 0.99999, FDR)) %>%
  inner_join(definitions) %>%
  filter(contrast %in% c("induction_only", "imipenem_partial")) %>%
  arrange(contrast, logFC) %>%
  group_by(contrast) %>%
  mutate(index_asc = row_number()) %>%
  arrange(contrast, desc(logFC)) %>%
  mutate(index_desc = row_number()) %>%
  ungroup() %>%
  mutate(gene_name = ifelse(is.na(gene_name), sprintf("bold('%s')", locus_tag), sprintf("italic('%s')", gene_name))) %>%
  mutate(type = ifelse(type %in% c("mismatch", "perfect essential"), "essential", "nonessential")) %>%
  ggplot(aes(y = logFC, x = FDR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "black") +
  scale_x_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
  geom_point(aes(fill = type, color = type), size = 2, alpha = 0.35, shape = 16) +
  scale_color_manual(
    values = c(
      # "control" = "#bbbbbb",
      "nonessential" = "#1F78B4",
      # "mismatch" = "#FF7F00",
      "essential" = "#E31A1C")) +
        scale_fill_manual(
    values = c(
      # "control" = "#bbbbbb",
      "nonessential" = "#1F78B4",
      # "mismatch" = "#FF7F00",
      "essential" = "#E31A1C")) +
  theme_bw() +
  # keep the y-axis the same for all plots
  facet_wrap(~contrast) +
  theme(strip.text.y = element_text(angle = 0)) +
  geom_text_repel(
    data = . %>%
      group_by(contrast) %>% 
      filter(FDR < 0.01 & (index_asc <= 10 | index_desc <= 10)),
      aes(label = gene_name),
    # size = 3,
    # segment.size = 0.2,
    segment.color = "black",
    segment.alpha = 0.5,
    box.padding = 0.5,
    point.padding = 0.5,
    force = 10,
    nudge_x = 0,
    nudge_y = 0,
    parse = TRUE,
  ) +
  geom_point(
    data = . %>%
      group_by(contrast) %>% 
      filter(FDR < 0.01 & (index_asc <= 10 | index_desc <= 10)),
    aes(color = type), 
    size = 4,
    shape = 21,
    alpha = 0.5,
    stroke = 1) +
    # x axis should say Confidence (FDR) and y axis should say Relatie Fitness Score (log2FC)
  labs(x = "Confidence (FDR)", y = "Relative Fitness Score (log2FC)", color = "Type", fill = "Type") 

print(volcano_plots)