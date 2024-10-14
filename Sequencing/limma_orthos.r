source("Sequencing/gamma_limma.r")

ortho_enrichments <- orthologs %>%
  inner_join(
    rbind(
      eco_enrichments,
      ecl_enrichments,
      kpn_enrichments
    )
  ) %>%
  select(HOG, genes_targeted, term, description) %>%
  group_by(genes_targeted, term, description) %>%
  ungroup() %>%
  select(HOG, term, description) %>%
  unique() %>%
  inner_join(orthologs) %>%
  select(term, description, assembly, locus_tag) %>%
  filter(!term %like% "CL:")

ortho_terms <- ortho_enrichments %>%
  select(term, description) %>%
  unique()

eco_ortho_spacers_for_terms <- ortho_enrichments %>%
  inner_join(
    eco_targets,
    na_matches = "never",
    relationship = "many-to-many"
  )

ecl_ortho_spacers_for_terms <- ortho_enrichments %>%
  inner_join(
    ecl_targets,
    na_matches = "never",
    relationship = "many-to-many"
  )

kpn_ortho_spacers_for_terms <- ortho_enrichments %>%
  inner_join(
    kpn_targets,
    na_matches = "never",
    relationship = "many-to-many"
  )

eco_ortho_spacers_in_sets <- split(
  eco_ortho_spacers_for_terms$spacer,
  eco_ortho_spacers_for_terms$term
)

ecl_ortho_spacers_in_sets <- split(
  ecl_ortho_spacers_for_terms$spacer,
  ecl_ortho_spacers_for_terms$term
)

kpn_ortho_spacers_in_sets <- split(
  kpn_ortho_spacers_for_terms$spacer,
  kpn_ortho_spacers_for_terms$term
)

# Find the index of the spacers in the DGEList object
eco_ortho_spacers_in_sets_index <- lapply(
  eco_ortho_spacers_in_sets,
  match,
  rownames(eco_dge)
)

ecl_ortho_spacers_in_sets_index <- lapply(
  ecl_ortho_spacers_in_sets,
  match,
  rownames(ecl_dge)
)

kpn_ortho_spacers_in_sets_index <- lapply(
  kpn_ortho_spacers_in_sets,
  match,
  rownames(kpn_dge)
)

# Perform competitive gene set enrichment analysis

## E. coli
eco_ortho_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = eco_v,
    index = eco_ortho_spacers_in_sets_index,
    design = eco_design_matrix,
    weights = eco_v$E %>% data.table(keep.rownames = "spacer") %>% select(spacer) %>% inner_join(eco_v_targets %>% select(spacer, weight) %>% unique()) %>% `$`(weight),
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = ortho_terms$term),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  left_join(ortho_terms)

## E. cloacae
ecl_ortho_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = ecl_v,
    index = ecl_ortho_spacers_in_sets_index,
    design = ecl_design_matrix,
    weights = ecl_v$E %>% data.table(keep.rownames = "spacer") %>% select(spacer) %>% inner_join(ecl_v_targets %>% select(spacer, weight) %>% unique()) %>% `$`(weight),
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = ortho_terms$term),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  inner_join(ortho_terms)

## K. pneumoniae
kpn_ortho_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    y = kpn_v,
    index = kpn_ortho_spacers_in_sets_index,
    design = kpn_design_matrix,
    weights = kpn_v$E %>% data.table(keep.rownames = "spacer") %>% select(spacer) %>% inner_join(kpn_v_targets %>% select(spacer, weight) %>% unique()) %>% `$`(weight),
    contrast = contrast_column
  ) |>
    data.table(keep.rownames = "term") |>
    mutate(
      term = factor(term, levels = ortho_terms$term),
      contrast = contrast_name
    )
  result
}) %>%
  do.call(rbind, .) |>
  rename(NGuides = NGenes) |>
  left_join(ortho_terms)


create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, ortho_enrichments, eco_ortho_sets[term == "GO:0032153", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ortho_enrichments, ecl_ortho_sets[term == "GO:0032153", ], 12, "E. cloacae")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, ortho_enrichments, kpn_ortho_sets[term == "GO:0032153", ], 12, "K. pneumoniae")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, ortho_enrichments, eco_ortho_sets[term == "GOCC:0030428", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ortho_enrichments, ecl_ortho_sets[term == "GOCC:0030428", ], 12, "E. cloacae")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, ortho_enrichments, kpn_ortho_sets[term == "GOCC:0030428", ], 12, "K. pneumoniae")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, ortho_enrichments, eco_ortho_sets[term == "GO:0045271", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ortho_enrichments, ecl_ortho_sets[term == "GO:0045271", ], 12, "E. cloacae")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, ortho_enrichments, kpn_ortho_sets[term == "GO:0045271", ], 12, "K. pneumoniae")

create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, ortho_enrichments, eco_ortho_sets[term == "GO:0030428", ], 12, "*E. coli*")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ortho_enrichments, ecl_ortho_sets[term == "GO:0030428", ], 12, "*E. cloacae*")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, ortho_enrichments, kpn_ortho_sets[term == "GO:0030428", ], 12, "*K. pneumoniae*")

# do GO:0043038
create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, ortho_enrichments, eco_ortho_sets[term == "GO:0043038", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ortho_enrichments, ecl_ortho_sets[term == "GO:0043038", ], 12, "E. cloacae")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, ortho_enrichments, kpn_ortho_sets[term == "GO:0043038", ], 12, "K. pneumoniae")


# do map03020
create_plot(eco_design_matrix, eco_full_with_zeros, eco_v_targets, ortho_enrichments, eco_ortho_sets[term == "map03020", ], 12, "E. coli")
create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ortho_enrichments, ecl_ortho_sets[term == "map03020", ], 12, "E. cloacae")
create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, ortho_enrichments, kpn_ortho_sets[term == "map03020", ], 12, "K. pneumoniae")

# create_plot(ecl_design_matrix, ecl_full_with_zeros, ecl_v_targets, ecl_enrichments, ecl_sets[term == "GO:0032153", ], 12, "Enterobacter")
# create_plot(kpn_design_matrix, kpn_full_with_zeros, kpn_v_targets, kpn_enrichments, kpn_sets[term == "GO:0032153", ], 12, "Klebsiella")
