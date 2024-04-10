require(pacman)

# Load the packages
p_load(
  plotly, data.table, scales, edgeR, statmod, poolr, ggtext, viridis, ggforce, igraph,
  pheatmap, svglite, ggplot2, ggrepel, RColorBrewer, tidyverse, magrittr, ggpubr, ggallin,
  Hmisc
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

drop_treated_uninduced <- function(df) {
  df |> filter((imipenem == 0 & induced == FALSE) | induced == TRUE)
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

ecl_1_counts <- load_counts("Sequencing/IMI_Ecloacae/ECLimiRDW1-12.tsv.gz") |> filter_artifacts()
ecl_2_counts <- load_counts("Sequencing/IMI_Ecloacae/ECL_2_imiRDW1-12.tsv.gz") |> filter_artifacts()

ecl_counts <- ecl_1_counts |>
  rbind(ecl_2_counts) |>
  group_by(sample, spacer) |>
  summarise(count = sum(count)) |>
  data.table()

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
  inner_join(eco_design)

# eco_full <- eco_full |>
#   drop_high_dose()


ecl_full <- ecl_counts |>
  inner_join(ecl_types) |>
  inner_join(ecl_design)

# ecl_full <- ecl_full |>
#   drop_high_dose()


kpn_full <- kpn_counts |>
  inner_join(kpn_types) |>
  inner_join(kpn_design)

# kpn_full <- kpn_full |>
#   drop_high_dose()

eco_full <- eco_full |>
  drop_treated_uninduced()

ecl_full <- ecl_full |>
  drop_treated_uninduced()

kpn_full <- kpn_full |>
  drop_treated_uninduced()
 

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
  # "imipenem_1x",
  # "imipenem_2x",
  "induced_imipenem_1x",
  "induced_imipenem_2x"
)

# create design matrix for edgeR with the groups
eco_design_matrix <- model.matrix(
  ~ factor(induced) + factor(imipenem),
  data = eco_names
) |>
  set_rownames(eco_names$verbose) |>
  set_colnames(contrast_levels)

ecl_design_matrix <- model.matrix(
  ~ factor(induced) + factor(imipenem),
  data = ecl_names
) |>
  set_rownames(ecl_names$verbose) |>
  set_colnames(contrast_levels)

kpn_design_matrix <- model.matrix(
  ~ factor(induced) + factor(imipenem),
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
  intercept = intercept,
  induced = induced,
  # imipenem_1x = imipenem_1x,
  # imipenem_2x = imipenem_2x,
  induced_imipenem_1x = induced_imipenem_1x,
  induced_imipenem_2x = induced_imipenem_2x,
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


create_volcano_plot <- function(results, title, print_plot = TRUE, label_top = 10) {
  volcano_plots <- results |>
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
      color = "black",
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
    ) +
    geom_text_repel(
      data = rbind(
        results |>
        arrange(type, ((logFC))) |>
        group_by(contrast, type) |>
        slice_min(
          n = label_top,
          order_by = ((logFC))
        ),
        results |>
        arrange(type, (logFC)) |>
        group_by(contrast, type) |>
        slice_max(
          n = label_top,
          order_by = (logFC)
          )
        
      ),
      aes(label = gene),
      size = 3,
      segment.color = "grey50",
      segment.size = 0.5,
      segment.alpha = 0.5
    )
  if (print_plot) {
    print(volcano_plots)
  } else {
    htmlwidgets::saveWidget(plotly::ggplotly(volcano_plots), paste0(title, "_volcano.html"))
  }
}

create_volcano_plot(eco_results, "E. coli", print_plot = TRUE)
create_volcano_plot(ecl_results, "E. cloacae", print_plot = TRUE)
create_volcano_plot(kpn_results, "K. pneumoniae", print_plot = TRUE)

eco_results_median <- eco_results |>
  group_by(contrast, type, locus_tag, gene) |>
  summarise(
    logFC = median(logFC),
    adj.P.Val = poolr::stouffer(adj.P.Val)$p
  )

ecl_results_median <- ecl_results |>
  group_by(contrast, type, locus_tag, gene) |>
  summarise(
    logFC = median(logFC),
    adj.P.Val = poolr::stouffer(adj.P.Val)$p
  )

kpn_results_median <- kpn_results |>
  group_by(contrast, type, locus_tag, gene) |>
  summarise(
    logFC = median(logFC),
    adj.P.Val = poolr::stouffer(adj.P.Val)$p
  )

create_volcano_plot(
  eco_results_median |> filter(
    contrast %in%
      c(
        "induced",
        "induced_imipenem_1x",
        "induced_imipenem_2x"
      )
  ),
  "E. coli (median)",
  print_plot = TRUE
)

create_volcano_plot(
  ecl_results_median |> filter(
    contrast %in%
      c(
        "induced",
        "induced_imipenem_1x",
        "induced_imipenem_2x"
      )
  ),
  "E. cloacae (median)",
  print_plot = TRUE
)

create_volcano_plot(
  kpn_results_median |> filter(
    contrast %in%
      c(
        "induced",
        "induced_imipenem_1x",
        "induced_imipenem_2x"
      )
  ),
  "K. pneumoniae (median)",
  print_plot = TRUE
)

create_volcano_plot(
  eco_results_median |>
    filter(type != "mismatch" & type != "control") |>
    filter(contrast %in% c("induced", "induced_imipenem_1x", "induced_imipenem_2x")) |>
    arrange(logFC), "E. coli (median)",
  print_plot = TRUE
)

create_volcano_plot(
  ecl_results_median |>
    filter(type != "mismatch" & type != "control") |>
    filter(contrast %in% c("induced", "induced_imipenem_1x", "induced_imipenem_2x")) |>
    arrange(logFC), "E. cloacae (median)",
  print_plot = TRUE
)

create_volcano_plot(
  kpn_results_median |>
    filter(type != "mismatch" & type != "control") |>
    filter(contrast %in% c("induced", "induced_imipenem_1x", "induced_imipenem_2x")) |>
    arrange(logFC), "K. pneumoniae (median)",
  print_plot = TRUE
)

# StringDB Operations
## Load StringDB

eco_string <- fread("Organisms/511145.protein.enrichment.terms.v12.0.txt.gz") |>
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) |>
  unique()

ecl_string <- fread("Organisms/STRG0A99CYX.protein.enrichment.terms.v12.0.txt.gz") |>
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) |>
  unique()

kpn_string <- fread("Organisms/STRG0A57EUR.protein.enrichment.terms.v12.0.txt.gz") |>
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) |>
  unique()

## Group by category, term, and description
create_gene_groups <- function(df) {
  gene_groups <- df |>
    group_by(category, term, description) |>
    summarise(
      gene_count = n(),
      locus_tag = list(sort(unique(locus_tag)))
    ) |>
    mutate(
      locus_tag_group = vapply(
        locus_tag,
        paste,
        collapse = ",",
        FUN.VALUE = character(1)
      )
    )
  return(gene_groups)
}

find_complete_terms <- function(gene_groups, targets) {
  gene_groups %>%
    unnest(locus_tag) %>%
    inner_join(targets %>% select(locus_tag) %>% unique()) %>%
    group_by(term, gene_count) %>%
    summarise(genes_targeted = n()) %>%
    filter(gene_count == genes_targeted)
}




eco_gene_groups <- create_gene_groups(eco_string)
ecl_gene_groups <- create_gene_groups(ecl_string)
kpn_gene_groups <- create_gene_groups(kpn_string)

eco_complete_terms <- find_complete_terms(eco_gene_groups, eco_targets)
ecl_complete_terms <- find_complete_terms(ecl_gene_groups, ecl_targets)
kpn_complete_terms <- find_complete_terms(kpn_gene_groups, kpn_targets)

# only perform enrichments where all genes are available
eco_gene_groups <- eco_gene_groups |>
  inner_join(eco_complete_terms)

ecl_gene_groups <- ecl_gene_groups |>
  inner_join(ecl_complete_terms)

kpn_gene_groups <- kpn_gene_groups |>
  inner_join(kpn_complete_terms)

# Enrichment groups mapping genes to terms

eco_enrichments <- eco_gene_groups |>
  ungroup() |>
  distinct(locus_tag_group, .keep_all = TRUE) |>
  unnest(locus_tag)

ecl_enrichments <- ecl_gene_groups |>
  ungroup() |>
  distinct(locus_tag_group, .keep_all = TRUE) |>
  unnest(locus_tag)

kpn_enrichments <- kpn_gene_groups |>
  ungroup() |>
  distinct(locus_tag_group, .keep_all = TRUE) |>
  unnest(locus_tag)

# Make a list of all the enrichments

## Make sure there are at least 4 guides in the enrichment set
eco_measurable_terms <- eco_enrichments |>
  inner_join(
    eco_targets,
    relationship = "many-to-many"
  ) |>
  unique() |>
  group_by(term) |>
  tally() |>
  filter(n >= 4) |>
  filter(
    !term %in% (eco_enrichments |>
      full_join(
        eco_targets,
        relationship = "many-to-many"
      ) |>
      filter(is.na(spacer)) |>
      pull(term) |>
      unique())
  ) |>
  unique() |>
  arrange(n) |>
  inner_join(
    eco_enrichments |>
      select(term, description) |>
      unique()
  )

ecl_measurable_terms <- ecl_enrichments |>
  inner_join(
    ecl_targets,
    relationship = "many-to-many"
  ) |>
  unique() |>
  group_by(term) |>
  tally() |>
  filter(n >= 4) |>
  filter(
    !term %in% (ecl_enrichments |>
      full_join(
        ecl_targets,
        relationship = "many-to-many"
      ) |>
      filter(is.na(spacer)) |>
      pull(term) |>
      unique())
  ) |>
  unique() |>
  arrange(n) |>
  inner_join(
    ecl_enrichments |>
      select(term, description) |>
      unique()
  )

kpn_measurable_terms <- kpn_enrichments |>
  inner_join(
    kpn_targets,
    relationship = "many-to-many"
  ) |>
  unique() |>
  group_by(term) |>
  tally() |>
  filter(n >= 4) |>
  filter(
    !term %in% (kpn_enrichments |>
      full_join(
        kpn_targets,
        relationship = "many-to-many"
      ) |>
      filter(is.na(spacer)) |>
      pull(term) |>
      unique())
  ) |>
  unique() |>
  arrange(n) |>
  inner_join(
    kpn_enrichments |>
      select(term, description) |>
      unique()
  )

## Unique terms
eco_unique_terms <- eco_measurable_terms |>
  select(term, description) |>
  unique()

ecl_unique_terms <- ecl_measurable_terms |>
  select(term, description) |>
  unique()

kpn_unique_terms <- kpn_measurable_terms |>
  select(term, description) |>
  unique()

## Spacers in these terms

eco_spacers_for_terms <- eco_measurable_terms |>
  inner_join(eco_enrichments, relationship = "many-to-many") |>
  inner_join(eco_targets, relationship = "many-to-many")

ecl_spacers_for_terms <- ecl_measurable_terms |>
  inner_join(ecl_enrichments, relationship = "many-to-many") |>
  inner_join(ecl_targets, relationship = "many-to-many")

kpn_spacers_for_terms <- kpn_measurable_terms |>
  inner_join(kpn_enrichments, relationship = "many-to-many") |>
  inner_join(kpn_targets, relationship = "many-to-many")

# Split the spacer column by term
eco_spacers_in_sets <- split(eco_spacers_for_terms$spacer, eco_spacers_for_terms$term)
ecl_spacers_in_sets <- split(ecl_spacers_for_terms$spacer, ecl_spacers_for_terms$term)
kpn_spacers_in_sets <- split(kpn_spacers_for_terms$spacer, kpn_spacers_for_terms$term)

# make control set
make_control_set <- function(targets_df) {
  control_df <- targets_df %>%
    filter(type == "control") %>%
    pull(spacer) %>%
    data.table(spacer = ., term = "control")

  split(control_df$spacer, control_df$term)
}

eco_control_spacers_in_sets <- make_control_set(eco_targets)
ecl_control_spacers_in_sets <- make_control_set(ecl_targets)
kpn_control_spacers_in_sets <- make_control_set(kpn_targets)


# Find the indices of each set of locus tags in rownames(dge)
eco_spacers_in_sets_index <- lapply(eco_spacers_in_sets, match, rownames(eco_dge))
ecl_spacers_in_sets_index <- lapply(ecl_spacers_in_sets, match, rownames(ecl_dge))
kpn_spacers_in_sets_index <- lapply(kpn_spacers_in_sets, match, rownames(kpn_dge))

# Find the indices of each control
eco_control_spacers_in_sets_index <- lapply(eco_control_spacers_in_sets, match, rownames(eco_dge))
ecl_control_spacers_in_sets_index <- lapply(ecl_control_spacers_in_sets, match, rownames(ecl_dge))
kpn_control_spacers_in_sets_index <- lapply(kpn_control_spacers_in_sets, match, rownames(kpn_dge))

# Voom again?
eco_v <- voomWithQualityWeights(eco_dge, eco_design_matrix, plot = TRUE)
ecl_v <- voomWithQualityWeights(ecl_dge, ecl_design_matrix, plot = TRUE)
kpn_v <- voomWithQualityWeights(kpn_dge, kpn_design_matrix, plot = TRUE)

eco_v_targets <- eco_v$E |>
  data.table(keep.rownames = "spacer") |>
  select(spacer) |>
  left_join(eco_targets)

ecl_v_targets <- ecl_v$E |>
  data.table(keep.rownames = "spacer") |>
  select(spacer) |>
  left_join(ecl_targets)

kpn_v_targets <- kpn_v$E |>
  data.table(keep.rownames = "spacer") |>
  select(spacer) |>
  left_join(kpn_targets)


# Weight based on y_pred seems to be working
eco_v_targets$y_pred <- as.numeric(eco_v_targets$y_pred)
eco_v_targets[is.na(target) | target == "None", weight := min(eco_v_targets$y_pred, na.rm = TRUE)]
eco_v_targets[spacer == target, weight := max(eco_v_targets$y_pred, na.rm = TRUE)]
eco_v_targets[type == "mismatch", weight := y_pred]
eco_v_targets$weight <- scales::rescale(as.numeric(eco_v_targets$weight), to = c(1, 100))

ecl_v_targets$y_pred <- as.numeric(ecl_v_targets$y_pred)
ecl_v_targets[is.na(target) | target == "None", weight := min(ecl_v_targets$y_pred, na.rm = TRUE)]
ecl_v_targets[spacer == target, weight := max(ecl_v_targets$y_pred, na.rm = TRUE)]
ecl_v_targets[type == "mismatch", weight := y_pred]
ecl_v_targets$weight <- scales::rescale(as.numeric(ecl_v_targets$weight), to = c(1, 100))

kpn_v_targets$y_pred <- as.numeric(kpn_v_targets$y_pred)
kpn_v_targets[is.na(target) | target == "None", weight := min(kpn_v_targets$y_pred, na.rm = TRUE)]
kpn_v_targets[spacer == target, weight := max(kpn_v_targets$y_pred, na.rm = TRUE)]
kpn_v_targets[type == "mismatch", weight := y_pred]
kpn_v_targets$weight <- scales::rescale(as.numeric(kpn_v_targets$weight), to = c(1, 100))

# Get sets for each organism
## E. coli
eco_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = eco_v,
    index = eco_spacers_in_sets_index,
    design = eco_design_matrix,
    weights = eco_v$E %>% data.table(keep.rownames = "spacer") %>% select(spacer) %>% inner_join(eco_v_targets %>% select(spacer, weight) %>% unique()) %>% `$`(weight),
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = eco_unique_terms$term),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  left_join(eco_unique_terms)

eco_control_set <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = eco_v,
    index = eco_control_spacers_in_sets_index,
    design = eco_design_matrix,
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = "control"),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  mutate(description = "Control")

## E. cloacae
ecl_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = ecl_v,
    index = ecl_spacers_in_sets_index,
    design = ecl_design_matrix,
    weights = ecl_v$E %>% data.table(keep.rownames = "spacer") %>% select(spacer) %>% inner_join(ecl_v_targets %>% select(spacer, weight) %>% unique()) %>% `$`(weight),
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = ecl_unique_terms$term),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  left_join(ecl_unique_terms)

ecl_control_set <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = ecl_v,
    index = ecl_control_spacers_in_sets_index,
    design = ecl_design_matrix,
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = "control"),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  mutate(description = "Control")

## K. pneumoniae
kpn_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = kpn_v,
    index = kpn_spacers_in_sets_index,
    design = kpn_design_matrix,
    weights = kpn_v$E %>% data.table(keep.rownames = "spacer") %>% select(spacer) %>% inner_join(kpn_v_targets %>% select(spacer, weight) %>% unique()) %>% `$`(weight),
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = kpn_unique_terms$term),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  left_join(kpn_unique_terms)

kpn_control_set <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = kpn_v,
    index = kpn_control_spacers_in_sets_index,
    design = kpn_design_matrix,
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = "control"),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  mutate(description = "Control")

# Add zero values to the full data
# create a full data frame with the verbose names and add zeros for missing spacers
eco_full_with_zeros <- eco_full_casted %>% 
  melt(id.vars = "spacer", variable.name = "verbose", value.name = "count") %>% 
  inner_join(eco_names) %>% 
  inner_join(eco_design) %>%
  inner_join(eco_targets %>% select(spacer) %>% unique())

ecl_full_with_zeros <- ecl_full_casted %>%
  melt(id.vars = "spacer", variable.name = "verbose", value.name = "count") %>%
  inner_join(ecl_names) %>%
  inner_join(ecl_design) %>%
  inner_join(ecl_targets %>% select(spacer) %>% unique())
  
kpn_full_with_zeros <- kpn_full_casted %>%
  melt(id.vars = "spacer", variable.name = "verbose", value.name = "count") %>%
  inner_join(kpn_names) %>%
  inner_join(kpn_design) %>%
  inner_join(kpn_targets %>% select(spacer) %>% unique())

## Create a plot for each organism
imipenem_color_map <- c(
  "0" = "grey50",
  "0.125" = "#33A02C",
  "0.25" = "#6A3D9A",
  "0.5" = "#FF7F00",
  "1" = "#E31A1C"
)

# Function to create a plot for each organism and each gene set enrichment
create_plot <- function(design_matrix, full_data, targets, enrichments, sets, limit, title) {
  set_name <- deparse(substitute(sets))

  plot_data <- full_data %>%
    # mutate(imipenem = ifelse(is.na(imipenem), "Stock", imipenem), induced = ifelse(is.na(induced), "Stock", induced)) %>%
    group_by(sample) %>%
    # mutate(imipenem = factor(imipenem, levels = c("Stock", "0", "0.125", "0.25", "0.5", "1"))) %>%
    # mutate(induced = factor(induced, levels = c("Stock", "FALSE", "TRUE"))) %>%
    mutate(cpm = cpm(count)) %>%
    ungroup() %>%
    inner_join(targets, by = "spacer", relationship = "many-to-many") %>%
    inner_join(enrichments, relationship = "many-to-many") %>%
    inner_join(
      rbind(
        sets %>%
          head(limit)
      )
    ) %>%
    group_by(factor(imipenem), factor(induced), term) %>%
    ungroup()

  group_evidence <- design_matrix %>%
    data.table(keep.rownames = "verbose") |>
    melt(
      variable.name = "contrast",
      value.name = "identity",
      id.vars = "verbose"
    ) %>%
    filter(identity == 1) %>%
    mutate(id = row_number()) %>%
    group_by(verbose) %>%
    filter(id == max(id)) %>%
    select(verbose, contrast)

  # get the number of locus tags in each set, not the number of guides, which is NGuides
  plot_data <- plot_data %>%
    group_by(factor(imipenem), factor(induced), term) %>%
    mutate(NGenes = locus_tag %>% unique() %>% length()) %>%
    ungroup() %>%
    inner_join(group_evidence)


  plot_data <- plot_data %>%
    mutate(all_genes = max(NGenes)) %>%
    mutate(all_guides = max(NGuides)) %>%
    mutate(description = stringr::str_wrap(description, width = 30)) %>%
    mutate(facet_title = paste0("**", term, "**", " — ", description, " — ", all_genes, " gene(s), ", all_guides, " guides" )) %>%
    # mutate(facet_title = sub("([^ \n]+)", "**\\1**", facet_title)) %>%
    mutate(facet_title = gsub("\n", "<br>", facet_title)) %>%
    mutate(facet_title = factor(facet_title, levels = facet_title %>% unique()))

  max_cpm <- max(plot_data$cpm)

  group_labels <- plot_data %>%
    arrange(facet_title) %>%
    mutate(group_title = paste(
      paste0(
        "**FDR** = ", signif(FDR, 2),
        ", *genes* = **",
        NGenes, "**(", NGuides, ")"
      ),
      paste0("**[", contrast, "]**"),
      sep = "<br>"
    )) %>%
    select(group_title, induced, imipenem, FDR, Direction, NGuides, NGenes, group) %>%
    unique() %>%
    mutate(cpm = max_cpm * 1.75) %>%
    mutate(mini_label = paste0(Direction, " (", signif(FDR, 2), ")")) %>%
    mutate(mini_label = ifelse(FDR <= 0.05, mini_label, "")) %>%
    mutate(mini_label = factor(mini_label, levels = mini_label %>% unique()))


  quantiles <- plot_data %>%
    group_by(imipenem, induced, facet_title, group) %>%
    arrange(facet_title) %>%
    summarise(cpm = wtd.quantile(cpm, weights = weight, probs = c(0.5)))

  print(quantiles)
  print(group_labels)



  ggplot(
    plot_data,
    aes(
      y = cpm,
      # x = factor(induced),
      x = group
    )
  ) +
    geom_tile(
      data = data.frame(induced = FALSE, group = "uninduced_0"),
      aes(x = 1, y = 0),
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
      alpha = 0.35,
      shape = 20,
    ) +
    geom_boxplot(
      aes(weight = weight),
      width = 0.9,
      outlier.shape = NA,
      position = position_dodge2(preserve = "single"),
      fill = NA
    ) +
    facet_wrap(~facet_title) +
    scale_size(range = c(0.25, 2.5)) +
    scale_fill_manual(values = imipenem_color_map) +
    scale_color_manual(values = imipenem_color_map) +
    labs(
      x = NULL,
      y = "Counts per Million",
      color = "Imipenem (µg/mL)",
      fill = "Imipenem (µg/mL)",
      size = "Predicted Weight"
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
      breaks = c(10^(0:5)),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    scale_x_discrete(labels = c("TRUE" = "Induced", "FALSE" = "Uninduced")) +
    theme_minimal() +
    theme(
      strip.text = element_markdown(),
      axis.text.x = element_text(size = rel(1.3), color = "black", angle = 45, hjust = 1),
      axis.title.y = element_text(size = rel(1.3), color = "black")
    ) +
    geom_text(
      data = group_labels,
      aes(label = mini_label),
      position = position_dodge2(width = 0.9, preserve = "single"),
      size = 2.5
    ) +
    geom_label(
      data = quantiles,
      aes(label = round(cpm, 2)),
      position = position_dodge2(width = 0.9, preserve = "single"),
      alpha = 0.75,
      size = 2.5
    ) +
    scale_x_discrete(labels = c("Uninduced 0x", "Induced 0x", "Induced 1x", "Induced 2x")) +

    ggtitle(title)
}



# Plot the controls for each organism
create_plot(
  eco_design_matrix,
  eco_full_with_zeros,
  eco_targets %>% filter(type == "control") %>% mutate(weight = 1),
  eco_targets %>% filter(type == "control") %>% mutate(term = "control", description = "Control"),
  eco_control_set %>% mutate(FDR = PValue),
  12,
  "E. coli"
)

create_plot(
  ecl_design_matrix,
  ecl_full_with_zeros,
  ecl_targets %>% filter(type == "control") %>% mutate(weight = 1),
  ecl_targets %>% filter(type == "control") %>% mutate(term = "control", description = "Control"),
  ecl_control_set %>% mutate(FDR = PValue),
  12,
  "Enterobacter"
)

create_plot(
  kpn_design_matrix,
  kpn_full_with_zeros,
  kpn_targets %>% filter(type == "control") %>% mutate(weight = 1),
  kpn_targets %>% filter(type == "control") %>% mutate(term = "control", description = "Control"),
  kpn_control_set %>% mutate(FDR = PValue),
  12,
  "Klebsiella"
)

# create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0015155", ], 12, "E. coli")
# create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0008509", ], 12, "E. coli")
# create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0008514", ], 12, "E. coli")
# create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "IPR014286", ], 12, "E. coli")
# create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0009245", ], 12, "E. coli")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0032153", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ecl_enrichments, ecl_sets[term == "GO:0032153", ], 12, "Enterobacter")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, kpn_enrichments, kpn_sets[term == "GO:0032153", ], 12, "Klebsiella")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0007231", ], 12, "E. coli")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0009254", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ecl_enrichments, ecl_sets[term == "GO:0009254", ], 12, "Enterobacter")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, kpn_enrichments, kpn_sets[term == "GO:0009254", ], 12, "Klebsiella")


create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GOCC:0030428", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ecl_enrichments, ecl_sets[term == "GOCC:0030428", ], 12, "Enterobacter")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, kpn_enrichments, kpn_sets[term == "GOCC:0030428", ], 12, "Klebsiella")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:1905153", ], 12, "E. coli")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GOCC:0045271", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ecl_enrichments, ecl_sets[term == "GOCC:0045271", ], 12, "Enterobacter")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, kpn_enrichments, kpn_sets[term == "GOCC:0045271", ], 12, "Klebsiella")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0032153", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ecl_enrichments, ecl_sets[term == "GO:0032153", ], 12, "Enterobacter")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, kpn_enrichments, kpn_sets[term == "GO:0032153", ], 12, "Klebsiella")

# draw plots for GO:0044391
create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0044391", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ecl_enrichments, ecl_sets[term == "GO:0044391", ], 12, "Enterobacter")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, kpn_enrichments, kpn_sets[term == "GO:0044391", ], 12, "Klebsiella")

# draw plots for IPR012338
create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "IPR012338", ], 12, "E. coli")
# does not exist in Enterobacter or Klebsiella. So find orthologs
create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, eco_enrichments, eco_sets[term == "GO:0010958", ], 12, "E. coli")


GO:0010958