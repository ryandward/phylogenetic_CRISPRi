require("pacman")

# Load the packages
p_load(
  plotly, data.table, scales, edgeR, statmod, poolr, ggtext, viridis, ggforce, igraph,
  pheatmap, svglite, ggplot2, ggrepel, RColorBrewer, tidyverse, magrittr, ggpubr, ggallin
)


data <- fread("Sequencing/IMI_Kpneumoniae/KPNimiRDW1-12.tsv.gz",
  header = FALSE, sep = "\t", col.names = c("sample", "spacer", "count")
) %>% filter(!(spacer %like% "\\*"))


# load the design
design <- fread(
  "Sequencing/IMI_Kpneumoniae/design.tsv",
  header = TRUE, sep = "\t"
)

data <- data %>% inner_join(design, by = "sample")


targets <- fread("Organisms/C_K_pneumoniae.targets.tsv.gz",
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
  mutate(
    cpm = 1e6 * count / sum(count)
  ) %>%
  data.table() %>%
  dcast(type + spacer + imipenem + induced ~ replicate, value.var = "cpm", fill = 0) %>%
  arrange(desc(type)) %>%
  ggplot(aes(x = A, y = B)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(
    aes(color = type),
    size = 1,
    alpha = 0.5,
    shape = 16
  ) +
  scale_color_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C"
    )
  ) +
  stat_cor(
    method = "pearson",
    aes(label = after_stat(label)),
    label.x = 0,
    label.y = 5,
    na.rm = TRUE
  ) +
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
  theme(
    strip.text.y = element_text(angle = 0)
  ) +
  ggtitle("Klebsiella pneumoniae")

print(replicate_plot)

# 0.5 looks bad, so we'll remove it
data <- data %>%
  filter(!(imipenem == 0.5))

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

colnames(design_matrix) <- c("intercept", "induced", "imipenem", "induced_imipenem")

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
  dcast(
    spacer ~ factor(
      verbose,
      levels = unique(verbose)
    ),
    value.var = "count"
  ) %>%
  select(-spacer) %>%
  data.matrix() %>%
  set_rownames(spacers)

# create a DGEList object
dge <- DGEList(
  counts = full_data_matrix,
  group = design_names$group,
  samples = design_names$verbose
)

# # normalize
dge <- calcNormFactors(dge)

# estimate dispersion
dge <- estimateGLMRobustDisp(dge, design_matrix)

# voom
v <- voomWithQualityWeights(dge, design_matrix)

# fit the model
voom_fit <- lmFit(v, design_matrix, method = "robust")

# these are the single-level "contrasts" we want to test
contrasts <- makeContrasts(
  intercept = intercept,
  induced = induced,
  induced_imipenem = induced_imipenem,
  imipenem = imipenem,
  levels = design_matrix
)

# create a single data.table with all the results, going through each level of the contrasts
results <- lapply(
  colnames(contrasts),
  function(level) {
    topTable(
      eBayes(
        contrasts.fit(
          voom_fit,
          contrasts = contrasts[, level]
        )
      ),
      n = Inf
    ) %>%
      data.table(keep.rownames = "spacer") %>%
      inner_join(targets) %>%
      mutate(contrast = level)
  }
) %>%
  rbindlist()

# factor the contrast
results$contrast <- factor(
  results$contrast,
  levels = colnames(contrasts)
)

# draw volcano plots with facets
volcano_plots <- results %>%
  filter(contrast %in% c("induced", "induced_imipenem")) %>%
  arrange(desc(type)) %>%
  ggplot(aes(y = logFC, x = adj.P.Val)) +
  scale_x_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
  geom_point(aes(color = type), size = 2, alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  scale_color_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C"
    )
  ) +
  theme_bw() +
  labs(color = "Type") +
  facet_wrap(contrast ~ ., scales = "free_x") +
  theme(strip.text.y = element_text(angle = 0))

print(volcano_plots)

################


all_string <- fread("Organisms/STRG0A57EUR.protein.enrichment.terms.v12.0.txt.gz") %>%
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
  unique()


gene_groups <- all_string %>%
  # filter(term %in% (all_string %>% group_by(term) %>% tally() %>% pull(unique(term)))) %>%
  group_by(category, term, description) %>%
  summarise(gene_count = n(), locus_tag = list(sort(unique(locus_tag)))) %>%
  mutate(locus_tag_group = vapply(locus_tag, paste, collapse = ",", FUN.VALUE = character(1)))

complete_terms <- gene_groups %>%
  unnest(locus_tag) %>%
  inner_join(
    targets %>%
      select(locus_tag) %>%
      unique()
  ) %>%
  group_by(term, gene_count) %>%
  summarize(genes_targeted = n()) %>%
  filter(gene_count == genes_targeted)

# only perform enrichments where all genes are available
gene_groups <- complete_terms %>% inner_join(gene_groups)

repeated_gene_groups <- gene_groups %>%
  group_by(locus_tag) %>%
  mutate(times_listed = n()) %>%
  arrange(locus_tag) %>%
  ungroup()

# pick the best annotation for each locus_tag_group, i.e., highest in term, and the lowest in the category_rank
ranked_annotations <- repeated_gene_groups %>%
  group_by(locus_tag_group, category) %>%
  arrange(versionsort::ver_sort(term)) %>%
  slice(n()) %>%
  ungroup() %>%
  mutate(category_rank = case_when(
    category == "Biological Process (Gene Ontology)" ~ 1,
    category == "Molecular Function (Gene Ontology)" ~ 2,
    category == "Cellular Component (Gene Ontology)" ~ 3,
    category == "Protein Domains and Features (InterPro)" ~ 4,
    category == "Protein Domains (SMART)" ~ 5,
    category == "Protein Domains (Pfam)" ~ 6,
    category == "Annotated Keywords (UniProt)" ~ 7,
    category == "Reactome Pathways" ~ 8,
    category == "Subcellular localization (COMPARTMENTS)" ~ 9,
    category == "Local Network Cluster (STRING)" ~ 10,
    TRUE ~ NA_integer_
  )) %>%
  group_by(locus_tag_group) %>%
  filter(category_rank == min(category_rank))


enrichments <- ranked_annotations %>%
  ungroup() %>%
  distinct(locus_tag_group, .keep_all = TRUE) %>%
  select(-locus_tag_group) %>%
  unnest(locus_tag)

enrichments_with_enough_spacers <- enrichments %>%
  inner_join(targets, relationship = "many-to-many") %>%
  unique() %>%
  group_by(term) %>%
  tally() %>%
  filter(n >= 4) %>%
  filter(!term %in% (enrichments %>% full_join(targets, relationship = "many-to-many") %>% filter(is.na(spacer)) %>% pull(term) %>% unique())) %>%
  unique() %>%
  arrange(n) %>%
  inner_join(enrichments %>% select(term, description) %>% unique())


# Get the unique terms
unique_terms <- unique(enrichments_with_enough_spacers$term)

target_spacers_for_terms <- enrichments_with_enough_spacers %>%
  inner_join(enrichments, relationship = "many-to-many") %>%
  inner_join(targets, relationship = "many-to-many")


#########################################################################################


# Split the spacer column by term
spacers_in_sets <- split(target_spacers_for_terms$spacer, target_spacers_for_terms$term)

# Find the indices of each set of locus tags in rownames(dge)
spacers_in_sets_index <- lapply(spacers_in_sets, function(locus_tags) which(rownames(dge) %in% locus_tags))

v <- voomWithQualityWeights(dge, design_matrix, plot = TRUE)

v_targets <- v$E %>%
  data.table(keep.rownames = "spacer") %>%
  select(spacer) %>%
  left_join(
    targets %>%
      filter(locus_tag %in% all_string$locus_tag) %>%
      group_by(spacer) %>%
      filter(
        is.na(target) |
          target == "None" |
          (sp_dir != tar_dir & abs(as.numeric(offset)) == min(abs(as.numeric(offset))) & overlap == max(overlap))
      ) %>%
      group_by(target)
  )

v_targets[y_pred == "None", y_pred := NA_integer_]
v_targets$y_pred <- as.numeric(v_targets$y_pred)
v_targets[is.na(target) | target == "None", weight := min(v_targets$y_pred, na.rm = TRUE)]
v_targets[spacer == target, weight := max(v_targets$y_pred, na.rm = TRUE)]
v_targets[type == "mismatch", weight := y_pred]
v_targets$weight <- rescale(as.numeric(v_targets$weight), to = c(1, 100))


# Perform the competitive gene set test for all gene sets
all_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v,
    index = spacers_in_sets_index,
    design = design_matrix,
    weights = v_targets$weight,
    # inter.gene.cor = 0.01,
    contrast = contrast_column
  ) %>%
    data.table(keep.rownames = "term") %>%
    mutate(term = factor(term, levels = unique_terms), contrast = contrast_name)
  result
}) %>%
  do.call(rbind, .) %>%
  rename(NGuides = NGenes)


############################################


create_plot <- function(full_data, targets, enrichments, sets, limit) {
  set_name <- deparse(substitute(sets))

  plot_data <- full_data %>%
    select(-type) %>%
    mutate(imipenem = ifelse(is.na(imipenem), "Stock", imipenem), induced = ifelse(is.na(induced), "Stock", induced)) %>%
    group_by(sample) %>%
    mutate(imipenem = factor(imipenem, levels = c("Stock", "0", "0.125", "0.25"))) %>%
    mutate(induced = factor(induced, levels = c("Stock", "FALSE", "TRUE"))) %>%
    mutate(cpm = cpm(count)) %>%
    ungroup() %>%
    inner_join(targets, by = "spacer", relationship = "many-to-many") %>%
    inner_join(enrichments, relationship = "many-to-many") %>%
    inner_join(
      rbind(
        sets %>%
          arrange(FDR) %>%
          head(limit)
      )
    ) %>%
    group_by(factor(imipenem), factor(induced), term) %>%
    ungroup() %>%
    arrange(FDR)

  # get the number of locus tags in each set, not the number of guides, which is NGuides
  plot_data <- plot_data %>%
    group_by(factor(imipenem), factor(induced), term) %>%
    mutate(NGenes = locus_tag %>% unique() %>% length()) %>%
    ungroup()

  plot_data <- plot_data %>%
    mutate(description = stringr::str_wrap(description, width = 30)) %>%
    mutate(
      facet_title = paste0("**", term, "**", " — ", Direction, "<br>", description)
    ) %>%
    # mutate(facet_title = sub("([^ \n]+)", "**\\1**", facet_title)) %>%
    mutate(
      facet_title = gsub("\n", "<br>", facet_title)
    ) %>%
    mutate(
      facet_title = paste(
        facet_title,
        paste0("**FDR** = ", signif(FDR, 2), ", *genes* = **", NGenes, "**(", NGuides, " guides)"),
        paste0("**[", contrast, "]**"),
        sep = "<br>"
      )
    ) %>%
    mutate(
      facet_title = factor(
        facet_title,
        levels = facet_title %>% unique()
      )
    )

  ggplot(plot_data, aes(y = cpm, x = factor(induced), group = interaction(factor(imipenem), factor(induced)))) +
    geom_tile(
      data = data.frame(induced = "TRUE"),
      aes(
        x = induced,
        y = 0
      ),
      width = 1,
      height = Inf,
      fill = "grey50",
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    geom_sina(
      aes(
        color = factor(imipenem),
        size = weight,
        weight = weight
      ),
      alpha = 0.5, shape = 16
    ) +
    geom_violin(
      aes(
        weight = as.numeric(weight)
      ),
      alpha = 0.25,
      draw_quantiles = c(0.25, 0.5, 0.75),
      position = "dodge"
    ) +
    facet_wrap(
      ~facet_title,
      nrow = 3,
      scales = "free_y"
    ) +
    scale_size(range = c(0.25, 2.5)) +
    # use some nice colors for fill and color
    scale_fill_manual(values = c("0" = "#1F78B4", "0.125" = "#FF7F00", "0.25" = "#E31A1C")) +
    scale_color_manual(values = c("0" = "#1F78B4", "0.125" = "#FF7F00", "0.25" = "#E31A1C")) +
    # bottom label should say induced and uninduced instead of true/false
    labs(x = NULL, y = "Counts per Million", color = "Imipenem (µg/mL)", fill = "Imipenem (µg/mL)", size = "Predicted Weight") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(10^(0:5)),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    scale_x_discrete(labels = c("TRUE" = "Induced", "FALSE" = "Uninduced")) +
    theme_minimal() +
    theme(
      strip.text = element_markdown(),
      axis.text.x = element_text(size = rel(1.3), color = "black", angle = 45, hjust = 1),
      axis.title.y = element_text(size = rel(1.3), color = "black")
    )
}

create_plot(full_data, v_targets, enrichments, all_sets[contrast == "induced_imipenem" & term == "CL:1801"], 12)
create_plot(full_data, v_targets, enrichments, all_sets[contrast == "induced_imipenem" & term == "CL:1801"], 12)
create_plot(full_data, v_targets, enrichments, all_sets[contrast == "induced_imipenem" & term == "GO:0006352"], 12)
create_plot(full_data, v_targets, enrichments, all_sets[contrast == "induced_imipenem" & term == "CL:1306"], 12)
create_plot(full_data, v_targets, enrichments, all_sets[contrast == "induced_imipenem" & term == "CL:1803"], 12)
create_plot(full_data, v_targets, enrichments, all_sets[contrast == "induced_imipenem" & term == "KW-0699"], 12)
