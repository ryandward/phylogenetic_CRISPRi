require(pacman)

# Load the packages
p_load(
  plotly, data.table, scales, edgeR, statmod, poolr, ggtext, viridis, ggforce, igraph,
  pheatmap, svglite, ggplot2, ggrepel, RColorBrewer, tidyverse, magrittr, ggpubr, ggallin
)

# Auxiliary functions
## fread_tsv
fread_tsv <- function(file, ...) {
  fread(
    file,
    header = TRUE,
    sep = "\t",
    na.strings = "None",
    ...
  )
}

## Load counts and  assign col.names 1: sample, 2: spacer, 3: count
load_counts <- function(file) {
  fread(
    file,
    header = FALSE,
    sep = "\t",
    col.names = c("sample", "spacer", "count"),
    na.strings = "None"
  )
}

## Filter out spacers with asterisks, which are potentially artifacts
filter_artifacts <- function(df) {
  df |> filter(!(spacer %like% "\\*"))
}

## Assign types based on mismatches
assign_types <- function(df) {
  df |>
    mutate(
      type = case_when(
        mismatches == 0 ~ "perfect",
        mismatches == 1 ~ "mismatch",
        is.na(mismatches) ~ "control"
      )
    )
}

## Standardize "target" column, mismatch is given in lower case, but we want to compare it with the reference genome so we need to convert it to upper case
standardize_targets <- function(df) {
  df |>
    mutate(
      target = toupper(target)
    )
}

## Mismatch guides were only generated for essential genes, so we can use them to identify essential genes
identify_essentials <- function(df) {
  essentials <- df |>
    filter(type == "mismatch") |>
    pull(locus_tag) |>
    unique()
}

## Extract types
extract_types <- function(df) {
  df |>
    select(spacer, chr, target, mismatches, type) |>
    unique()
}

## Drop high dose
drop_high_dose <- function(df) {
  df |> filter(imipenem != max(imipenem))
}

# Auxiliary mappings
color_type_map <- c(
  "control" = "#bbbbbb",
  "perfect" = "#1F78B4",
  "mismatch" = "#FF7F00",
  "perfect essential" = "#E31A1C"
)

# Load N0.tsv orthologs
orthologs <- fread_tsv("Sequencing/Orthology_N0.tsv") |>
  melt(
    id.vars = c("HOG", "OG", "Gene Tree Parent Clade"),
    variable.name = "assembly",
    value.name = "protein"
  ) |>
  filter(!is.na(protein) & protein != "") |>
  separate_wider_delim(
    cols = protein,
    delim = " ",
    names = c("locus_tag", "protein"),
    too_many = "merge",
    too_few = "align_start"
  )

# Load the data, first column is sample, second column is spacer, third column is count
eco_counts <- load_counts("Sequencing/IMI_Ecoli/ECimiRDW1-12.tsv.gz") |> filter_artifacts()
ecl_counts <- load_counts("Sequencing/IMI_Ecloacae/ECLimiRDW1-12.tsv.gz") |> filter_artifacts()
kpn_counts <- load_counts("Sequencing/IMI_Kpneumoniae/KPNimiRDW1-12.tsv.gz") |> filter_artifacts()

# Load experimental design
eco_design <- fread_tsv("Sequencing/IMI_Ecoli/design.tsv")
ecl_design <- fread_tsv("Sequencing/IMI_Ecloacae/design.tsv")
kpn_design <- fread_tsv("Sequencing/IMI_Kpneumoniae/design.tsv")

# Load target annotations
eco_targets <- fread_tsv("Organisms/A_E_coli.targets.tsv.gz") |>
  assign_types() |>
  standardize_targets()

ecl_targets <- fread_tsv("Organisms/B_E_cloacae.targets.tsv.gz") |>
  assign_types() |>
  standardize_targets()

kpn_targets <- fread_tsv("Organisms/C_K_pneumoniae.targets.tsv.gz") |>
  assign_types() |>
  standardize_targets()

# Identify essential genes
eco_essentials <- eco_targets |> identify_essentials()
ecl_essentials <- ecl_targets |> identify_essentials()
kpn_essentials <- kpn_targets |> identify_essentials()

# Provide a more nuanced classification based on essentiality
# Use data.table to update by reference

### E. coli

eco_targets[locus_tag %in% eco_essentials & type == "perfect", type := "perfect essential"]
eco_targets[spacer %in% eco_targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]
eco_types <- eco_targets |> extract_types()

### E. cloacae

ecl_targets[locus_tag %in% ecl_essentials & type == "perfect", type := "perfect essential"]
ecl_targets[spacer %in% ecl_targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]
ecl_types <- ecl_targets |> extract_types()

### K. pneumoniae

kpn_targets[locus_tag %in% kpn_essentials & type == "perfect", type := "perfect essential"]
kpn_targets[spacer %in% kpn_targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]
kpn_types <- kpn_targets |> extract_types()

# Replicate plots

generate_replicate_plot <- function(data, title) {
  cpm_data <- data |>
    group_by(sample) |>
    mutate(cpm = 1e6 * count / sum(count))

  relative_type_abundance <- cpm_data |>
    group_by(type) |>
    tally() |>
    arrange(desc(n))

  max_cpm <- cpm_data |>
    filter(!is.na(cpm)) |>
    pull(cpm) |>
    max()

  cpm_block <- cpm_data |>
    data.table() |>
    dcast(type + spacer + imipenem + induced ~ replicate, value.var = "cpm", fill = 0) |>
    inner_join(relative_type_abundance, by = "type") |>
    arrange(desc(n)) |>
    mutate(
      type = fct_reorder(type, n)
    )

  replicate_plot <- cpm_block |>
    ggplot(aes(x = A, y = B)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50",
      lwd = 1
    ) +
    geom_point(
      aes(color = type),
      size = 1,
      alpha = 0.5,
      shape = 16
    ) +
    scale_color_manual(values = color_type_map) +
    stat_cor(
      method = "pearson",
      aes(label = after_stat(label)),
      label.x = 0.5,
      label.y = log10(max_cpm) - 0.5,
      na.rm = TRUE
    ) +
    scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(10^(0:6)),
      labels = scales::label_number(scale_cut = scales::cut_long_scale()),
      limits = c(0, max_cpm)
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(10^(0:6)),
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      limits = c(0, max_cpm)
    ) +
    facet_grid(
      induced ~ imipenem,
      labeller = labeller(
        .cols = function(rows) paste("Imipenem:", rows),
        .rows = function(cols) ifelse(cols == "TRUE", "Induced", "Not induced")
      )
    ) +
    theme_bw() +
    guides(
      color = guide_legend(override.aes = list(size = 4))
    ) +
    labs(
      color = "**Type**",
      title = title
    ) +
    theme(
      legend.title = element_markdown(),
      plot.title = element_markdown(),
      strip.text.y = element_text(angle = 0),
      strip.text = element_markdown()
    )
  print(replicate_plot)
}

generate_replicate_plot(
  eco_counts |>
    inner_join(eco_types) |>
    inner_join(eco_design),
  "*E. coli* [**0x**, **1x**, **2x** imipenem]"
)

generate_replicate_plot(
  ecl_counts |>
    inner_join(ecl_types) |>
    inner_join(ecl_design),
  "*Enterobacter cloacae* [**0x**, **1x**, **2x** imipenem]"
)

generate_replicate_plot(
  kpn_counts |>
    inner_join(kpn_types) |>
    inner_join(kpn_design),
  "*Klebsiella pneumoniae* [**0x**, **1x**, **2x** imipenem]"
)
# Annotate the count data with experimental design and type
eco_full <- eco_counts |>
  inner_join(eco_types) |>
  inner_join(eco_design) |>
  drop_high_dose()

ecl_full <- ecl_counts |>
  inner_join(ecl_types) |>
  inner_join(ecl_design) |>
  drop_high_dose()

kpn_full <- kpn_counts |>
  inner_join(kpn_types) |>
  inner_join(kpn_design) |>
  drop_high_dose()

# Create names for the samples
eco_names <- eco_full |>
  select(sample, induced, imipenem, replicate) |>
  unique() |>
  mutate(
    induction = ifelse(induced == "TRUE", "induced", "uninduced")
  ) |>
  arrange(imipenem, desc(induction), replicate) |>
  mutate(
    verbose = paste(induction, imipenem, replicate),
    group = paste(induction, imipenem, sep = "_")
  ) |>
  mutate(group = gsub("\\.", "_", group)) |>
  select(-induction)

ecl_names <- ecl_full |>
  select(sample, induced, imipenem, replicate) |>
  unique() |>
  mutate(
    induction = ifelse(induced == "TRUE", "induced", "uninduced")
  ) |>
  arrange(imipenem, desc(induction), replicate) |>
  mutate(
    verbose = paste(induction, imipenem, replicate),
    group = paste(induction, imipenem, sep = "_")
  ) |>
  mutate(group = gsub("\\.", "_", group)) |>
  select(-induction)

kpn_names <- kpn_full |>
  select(sample, induced, imipenem, replicate) |>
  unique() |>
  mutate(
    induction = ifelse(induced == "TRUE", "induced", "uninduced")
  ) |>
  arrange(imipenem, desc(induction), replicate) |>
  mutate(
    verbose = paste(induction, imipenem, replicate),
    group = paste(induction, imipenem, sep = "_")
  ) |>
  mutate(group = gsub("\\.", "_", group)) |>
  select(-induction)

# create factor for the design
eco_names$group <- factor(eco_names$group, levels = unique(eco_names$group))
ecl_names$group <- factor(ecl_names$group, levels = unique(ecl_names$group))
kpn_names$group <- factor(kpn_names$group, levels = unique(kpn_names$group))

# explicitly reorder design_names by the groups (probably not necessary, but just in case)
eco_names <- eco_names |> arrange(group)
ecl_names <- ecl_names |> arrange(group)
kpn_names <- kpn_names |> arrange(group)

# then factor the design_names$sample by existing order, i.e., the order of the design_names$group
eco_names$sample <- factor(eco_names$sample, levels = unique(eco_names$sample))
ecl_names$sample <- factor(ecl_names$sample, levels = unique(ecl_names$sample))
kpn_names$sample <- factor(kpn_names$sample, levels = unique(kpn_names$sample))

# Experiment Level Names
contrast_levels <- c(
  "intercept",
  "induced",
  "imipenem_x1",
  "induced_imipenem_x1"
)

# create design matrix for edgeR with the groups
eco_design_matrix <- model.matrix(
  ~ factor(induced) * factor(imipenem),
  data = eco_names
) |>
  set_rownames(eco_names$verbose) |>
  set_colnames(contrast_levels)

ecl_design_matrix <- model.matrix(
  ~ factor(induced) * factor(imipenem),
  data = ecl_names
) |>
  set_rownames(ecl_names$verbose) |>
  set_colnames(contrast_levels)

kpn_design_matrix <- model.matrix(
  ~ factor(induced) * factor(imipenem),
  data = kpn_names
) |>
  set_rownames(kpn_names$verbose) |>
  set_colnames(contrast_levels)

# explicitly reorder the full data by the design groups
eco_full <- eco_names |>
  inner_join(eco_full) |>
  arrange(group)

ecl_full <- ecl_names |>
  inner_join(ecl_full) |>
  arrange(group)

kpn_full <- kpn_names |>
  inner_join(kpn_full) |>
  arrange(group)

# list of spacer names
# eco_spacers <- eco_full |>
#   pull(spacer) |>
#   unique()
# ecl_spacers <- ecl_full |>
#   pull(spacer) |>
#   unique()
# kpn_spacers <- kpn_full |>
#   pull(spacer) |>
#   unique()

# create a matrix of counts for edgeR
## E. coli
eco_full_casted <- eco_full |>
  dcast(
    spacer ~ factor(
      verbose,
      levels = unique(verbose)
    ),
    value.var = "count", fill = 0
  )

eco_full_matrix <- eco_full_casted |>
  select(-spacer) |>
  data.matrix() |>
  set_rownames(eco_full_casted$spacer)


## E. cloacae
ecl_full_casted <- ecl_full |>
  dcast(
    spacer ~ factor(
      verbose,
      levels = unique(verbose)
    ),
    value.var = "count", fill = 0
  )

ecl_full_matrix <- ecl_full_casted |>
  select(-spacer) |>
  data.matrix() |>
  set_rownames(ecl_full_casted$spacer)

## K. pneumoniae
kpn_full_casted <- kpn_full |>
  dcast(
    spacer ~ factor(
      verbose,
      levels = unique(verbose)
    ),
    value.var = "count", fill = 0
  )

kpn_full_matrix <- kpn_full_casted |>
  select(-spacer) |>
  data.matrix() |>
  set_rownames(kpn_full_casted$spacer)

# create DGEList object
eco_dge <- DGEList(counts = eco_full_matrix, group = eco_names$group)
ecl_dge <- DGEList(counts = ecl_full_matrix, group = ecl_names$group)
kpn_dge <- DGEList(counts = kpn_full_matrix, group = kpn_names$group)

# normalize the DGEList object
eco_dge <- calcNormFactors(eco_dge)
ecl_dge <- calcNormFactors(ecl_dge)
kpn_dge <- calcNormFactors(kpn_dge)

# estimate dispersion
eco_dge <- estimateDisp(eco_dge, design = eco_design_matrix)
ecl_dge <- estimateDisp(ecl_dge, design = ecl_design_matrix)
kpn_dge <- estimateDisp(kpn_dge, design = kpn_design_matrix)

# voom
eco_voom <- voom(eco_dge, design = eco_design_matrix)
ecl_voom <- voom(ecl_dge, design = ecl_design_matrix)
kpn_voom <- voom(kpn_dge, design = kpn_design_matrix)

# Fit the linear model
eco_fit <- lmFit(eco_voom, design = eco_design_matrix, method = "robust")
ecl_fit <- lmFit(ecl_voom, design = ecl_design_matrix, method = "robust")
kpn_fit <- lmFit(kpn_voom, design = kpn_design_matrix, method = "robust")

# common contrasts for all organisms
contrasts <- makeContrasts(
  "induced" = "induced",
  "imipenem_x1" = "imipenem_x1",
  "induced_imipenem_x1" = "induced_imipenem_x1",
  levels = contrast_levels
)

eco_results <- lapply(colnames(contrasts), function(contrast) {
  fit <- contrasts.fit(eco_fit, contrasts = contrasts[, contrast])
  topTable(eBayes(fit), n = Inf) |>
    data.table(keep.rownames = "spacer") |>
    inner_join(eco_targets) |>
    mutate(contrast = contrast)
}) |>
  rbindlist()

ecl_results <- lapply(colnames(contrasts), function(contrast) {
  fit <- contrasts.fit(ecl_fit, contrasts = contrasts[, contrast])
  topTable(eBayes(fit), n = Inf) |>
    data.table(keep.rownames = "spacer") |>
    inner_join(ecl_targets) |>
    mutate(contrast = contrast)
}) |>
  rbindlist()

kpn_results <- lapply(colnames(contrasts), function(contrast) {
  fit <- contrasts.fit(kpn_fit, contrasts = contrasts[, contrast])
  topTable(eBayes(fit), n = Inf) |>
    data.table(keep.rownames = "spacer") |>
    inner_join(kpn_targets) |>
    mutate(contrast = contrast)
}) |>
  rbindlist()



create_volcano_plot <- function(results, title, print_plot = TRUE) {
  volcano_plots <- results |>
    filter(
      contrast %in% c(
        "induced",
        "induced_imipenem_x1"
      )
    ) |>
    arrange(desc(type)) |>
    ggplot(
      aes(
        x = logFC,
        y = adj.P.Val,
        text = gene
      )
    ) +
    scale_y_continuous(
      trans = scales::reverse_trans() %of% scales::log10_trans()
    ) +
    geom_point(
      aes(color = type),
      size = 2,
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      color = "black",
      lwd = 0.5
    ) +
    geom_vline(
      xintercept = -1,
      linetype = "dashed",
      color = "#814b4b",
      lwd = 0.5
    ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "black",
      lwd = 0.5
    ) +
    scale_color_manual(
      values = color_type_map
    ) +
    facet_wrap(
      contrast ~ .,
      scales = "free_x"
    ) +
    theme_bw() +
    guides(
      color = guide_legend(override.aes = list(size = 4))
    ) +
    labs(
      color = "**Type**",
      title = title
    ) +
    theme(
      legend.title = element_markdown(),
      plot.title = element_markdown(),
      strip.text.y = element_text(angle = 0),
      strip.text = element_markdown()
    )
  if (print_plot) {
    print(volcano_plots)
  } else {
    htmlwidgets::saveWidget(plotly::ggplotly(volcano_plots), paste0(title, "_volcano.html"))
  }
}

create_volcano_plot(eco_results, "E. coli", print_plot = FALSE)
create_volcano_plot(ecl_results, "E. cloacae", print_plot = FALSE)
create_volcano_plot(kpn_results, "K. pneumoniae", print_plot = FALSE)
