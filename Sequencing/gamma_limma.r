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
  df %>% filter(!(spacer %like% "\\*"))
}

## Assign types based on mismatches
assign_types <- function(df) {
  df %>%
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
  df %>%
    mutate(
      target = toupper(target)
    )
}

## Mismatch guides were only generated for essential genes, so we can use them to identify essential genes
identify_essentials <- function(df) {
  essentials <- df %>%
    filter(type == "mismatch") %>%
    pull(locus_tag) %>%
    unique()

  df[locus_tag %in% essentials & type == "perfect", type := "perfect essential"]
  df[spacer %in% df[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]
  df
}

## Extract types
extract_types <- function(df) {
  df %>%
    select(spacer, chr, target, mismatches, type) %>%
    unique()
}

# Load N0.tsv orthologs
orthologs <- fread_tsv("Sequencing/Orthology_N0.tsv") %>%
  melt(
    id.vars = c("HOG", "OG", "Gene Tree Parent Clade"),
    variable.name = "assembly",
    value.name = "protein"
  ) %>%
  filter(!is.na(protein) & protein != "") %>%
  separate_wider_delim(
    cols = protein,
    delim = " ",
    names = c("locus_tag", "gene", "protein"),
    too_many = "merge",
    too_few = "align_start"
  )

# Load the data, first column is sample, second column is spacer, third column is count
eco_counts <- load_counts("Sequencing/IMI_Ecoli/ECimiRDW1-12.tsv.gz") %>% filter_artifacts()
ecl_counts <- load_counts("Sequencing/IMI_Ecloacae/ECLimiRDW1-12.tsv.gz") %>% filter_artifacts()
kpn_counts <- load_counts("Sequencing/IMI_Kpneumoniae/KPNimiRDW1-12.tsv.gz") %>% filter_artifacts()

# Load experimental design
eco_design <- fread_tsv("Sequencing/IMI_Ecoli/design.tsv")
ecl_design <- fread_tsv("Sequencing/IMI_Ecloacae/design.tsv")
kpn_design <- fread_tsv("Sequencing/IMI_Kpneumoniae/design.tsv")

# Load target annotations
eco_targets <- fread_tsv("Organisms/A_E_coli.targets.tsv.gz") %>%
  assign_types() %>%
  standardize_targets()

ecl_targets <- fread_tsv("Organisms/B_E_cloacae.targets.tsv.gz") %>%
  assign_types() %>%
  standardize_targets()

kpn_targets <- fread_tsv("Organisms/C_K_pneumoniae.targets.tsv.gz") %>%
  assign_types() %>%
  standardize_targets()

# Identify essential genes
eco_essentials <- eco_targets %>% identify_essentials()
ecl_essentials <- ecl_targets %>% identify_essentials()
kpn_essentials <- kpn_targets %>% identify_essentials()

# Provide a more nuanced classification based on essentiality
# Use data.table to update by reference

### E. coli

eco_targets[locus_tag %in% eco_essentials & type == "perfect", type := "perfect essential"]
eco_targets[spacer %in% eco_targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]
eco_types <- eco_targets %>% extract_types()

### E. cloacae

ecl_targets[locus_tag %in% ecl_essentials & type == "perfect", type := "perfect essential"]
ecl_targets[spacer %in% ecl_targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]
ecl_types <- ecl_targets %>% extract_types()

### K. pneumoniae

kpn_targets[locus_tag %in% kpn_essentials & type == "perfect", type := "perfect essential"]
kpn_targets[spacer %in% kpn_targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]
kpn_types <- kpn_targets %>% extract_types()

# Replicate plots

generate_replicate_plot <- function(data, title) {
  cpm_data <- data %>%
    group_by(sample) %>%
    mutate(cpm = 1e6 * count / sum(count))

  max_cpm <- cpm_data %>%
    filter(!is.na(cpm)) %>%
    pull(cpm) %>%
    max()

  cpm_block <- cpm_data %>%
    data.table() %>%
    dcast(type + spacer + imipenem + induced ~ sample, value.var = "cpm", fill = 0) %>%
    arrange(desc(type))

  replicate_plot <- cpm_block %>%
    ggplot(aes(x = A, y = B)) + # Ensure these aes mappings match your dataset.
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
    stat_cor(method = "pearson", aes(label = after_stat(label)), label.x = 0, label.y = log10(max_cpm) - 0.5, na.rm = TRUE) +
    scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(1, 10^(1:6)),
      labels = scales::label_number(scale_cut = scales::cut_long_scale()),
      limits = c(0, max_cpm)
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(1, 10^(1:6)),
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      limits = c(0, max_cpm)
    ) +
    facet_grid(imipenem ~ induced, labeller = labeller(
      .rows = function(rows) paste("Imipenem:", rows),
      .cols = function(cols) ifelse(cols == "TRUE", "Induced", "Not induced")
    )) +
    theme_bw() +
    labs(color = "Type") +
    theme(
      strip.text.y = element_text(angle = 0),
      strip.text = element_markdown()
    ) +
    labs(title = title) +
    theme(plot.title = element_markdown())

  print(replicate_plot)
}

generate_replicate_plot(
  eco_counts %>%
    inner_join(eco_types) %>%
    inner_join(eco_design), "*E. coli* [**0x**, **1x**, **2x** imipenem]"
)

generate_replicate_plot(
  ecl_counts %>%
    inner_join(ecl_types) %>%
    inner_join(ecl_design), "*Enterobacter cloacae* [**0x**, **1x**, **2x** imipenem]"
)

generate_replicate_plot(
  kpn_counts %>%
    inner_join(kpn_types) %>%
    inner_join(kpn_design), "*Klebsiella pneumoniae* [**0x**, **1x**, **2x** imipenem]"
)


