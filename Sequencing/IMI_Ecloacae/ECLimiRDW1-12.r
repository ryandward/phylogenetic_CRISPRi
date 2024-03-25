require("pacman")

# Load the packages
p_load(
  plotly, data.table, scales, edgeR, statmod, poolr, ggtext, viridis, ggforce, igraph,
  pheatmap, svglite, ggplot2, ggrepel, RColorBrewer, tidyverse, magrittr, ggpubr, ggallin
)


data <- fread("Sequencing/IMI_Ecloacae/ECLimiRDW1-12.tsv.gz",
  header = FALSE, sep = "\t", col.names = c("sample", "spacer", "count")
) %>% filter(!(spacer %like% "\\*"))


# load the design
design <- fread(
  "Sequencing/IMI_Ecloacae/design.tsv",
  header = TRUE, sep = "\t"
)

data <- data %>% inner_join(design, by = "sample")


targets <- fread("Organisms/B_E_cloacae.targets.tsv.gz",
  header = TRUE, sep = "\t", na.strings = "None"
) %>%
  mutate(
    type = case_when(
      mismatches == 0 ~ "perfect",
      mismatches == 1 ~ "mismatch",
      is.na(mismatches) ~ "control"
    )
  )

targets[type != "control", target := toupper(target)]


essentials <- targets %>%
  filter(type == "mismatch") %>%
  pull(locus_tag) %>%
  unique()

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
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = type), size = 1, alpha = 0.5, shape = 16) +
  scale_color_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C"
    )
  ) +
  stat_cor(method = "pearson", aes(label = after_stat(label)), label.x = 0, label.y = 5, na.rm = TRUE) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
    breaks = c(0, 1, 10^(1:6))
  ) +
  # labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
    breaks = c(0, 1, 10^(1:6))
  ) +
  # labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  facet_grid(imipenem ~ induced, labeller = labeller(
    .rows = function(rows) paste("Imipenem:", rows),
    .cols = function(cols) ifelse(cols == "TRUE", "Induced", "Not induced")
  )) +
  theme_bw() +
  labs(color = "Type") +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle("Enterobacter cloacae")

print(replicate_plot)

# 1.00 looks bad, so we'll remove it
data <- data %>% filter(!(imipenem == 1))

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
  select(sample, induced, imipenem, replicate) %>%
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
design_matrix <- model.matrix(~ factor(induced) * factor(imipenem), data = design_names) %>%
  set_rownames(design_names$verbose)

# colnames(design_matrix) <- c("intercept", "induced", "imipenem", "induced_imipenem")


# # explicitly reorder the full_data by the design_groups
# full_data <- design_names %>%
#   inner_join(full_data) %>%
#   arrange(group)


# # full_data <- full_data %>% inner_join(median_spacers)

# # list of spacer names
# spacers <- full_data %>%
#   select(spacer) %>%
#   unique() %>%
#   pull()

# # create a matrix of counts for edgeR
# full_data_matrix <- full_data %>%
#   dcast(spacer ~ factor(verbose, levels = unique(verbose)), value.var = "count") %>%
#   select(-spacer) %>%
#   data.matrix() %>%
#   set_rownames(spacers)

# # create a DGEList object
# dge <- DGEList(counts = full_data_matrix, group = design_names$group, samples = design_names$verbose)

# # # normalize
# dge <- calcNormFactors(dge)

# # estimate dispersion
# dge <- estimateGLMRobustDisp(dge, design_matrix)

# # voom
# v <- voomWithQualityWeights(dge, design_matrix)

# # fit the model
# voom_fit <- lmFit(v, design_matrix, method = "robust")

# # voom_fit <- voomLmFit(dge, design_matrix, sample.weights = TRUE, plot = TRUE)

# # create a list of contrasts
# contrasts <- makeContrasts(
#   # intercept = intercept,
#   induced = induced,
#   imipenem = imipenem,
#   induced_imipenem = induced_imipenem - imipenem,
#   imipenem_isolated = induced_imipenem - induced - imipenem,
#   levels = design_matrix
# )

# # create a single data.table with all the results, going through each contrast, one at a time
# results <- lapply(colnames(contrasts), function(contrast) {
#   fit <- contrasts.fit(voom_fit, contrasts = contrasts[, contrast])

#   topTable(eBayes(fit), n = Inf) %>%
#     # use_series(table) %>%
#     data.table(keep.rownames = "spacer") %>%
#     inner_join(targets) %>%
#     mutate(contrast = contrast)
# }) %>%
#   rbindlist()

# # add the coef as a factor, so that the order is preserved
# results$contrast <- factor(results$contrast, levels = colnames(contrasts))

# # draw volcano plots with facets
# volcano_plots <- results %>%
#   # filter(contrast %in% c("induced", "induced_imipenem")) %>%
#   arrange(desc(type)) %>%
#   ggplot(aes(x = logFC, y = adj.P.Val)) +
#   scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
#   geom_point(aes(color = type), size = 2, alpha = 0.5) +
#   # draw lines at 0.05  adj.P.Val and abs(logFC) of 1
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", size = 0.5) +
#   geom_vline(xintercept = -1, linetype = "dashed", color = "#814b4b", size = 0.5) +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
#   scale_color_manual(
#     values = c(
#       "control" = "#bbbbbb",
#       "perfect" = "#1F78B4",
#       "mismatch" = "#FF7F00",
#       "perfect essential" = "#E31A1C"
#     )
#   ) +
#   theme_bw() +
#   labs(color = "Type") +
#   facet_wrap(contrast ~ ., scales = "free_x") +
#   theme(strip.text.y = element_text(angle = 0))

# print(volcano_plots)


colnames(design_matrix) <- c("intercept", "induced", "imipenem", "induced_imipenem")

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
  dcast(spacer ~ factor(verbose, levels = unique(verbose)), value.var = "count") %>%
  select(-spacer) %>%
  data.matrix() %>%
  set_rownames(spacers)

# create a DGEList object
dge <- DGEList(counts = full_data_matrix, group = design_names$group, samples = design_names$verbose)

# # normalize
dge <- calcNormFactors(dge)

# estimate dispersion
dge <- estimateGLMRobustDisp(dge, design_matrix)

# fit the glmQLFTest
fit <- glmQLFit(dge, design = design_matrix)

contrasts <- makeContrasts(
  induced = induced,
  induced_imipenem = induced_imipenem,
  imipenem = imipenem,
  imipenem_difference = induced_imipenem - induced,
  # imipenem_extra = induced_imipenem - induced,
  levels = design_matrix
)

# create a single data.table with all the results, going through each contrast, one at a time
results <- lapply(colnames(contrasts), function(contrast) {
  topTags(glmQLFTest(fit, contrast = contrasts[, contrast]), n = Inf) %>%
    use_series(table) %>%
    data.table(keep.rownames = "spacer") %>%
    inner_join(targets) %>%
    mutate(contrast = contrast)
}) %>%
  rbindlist()

