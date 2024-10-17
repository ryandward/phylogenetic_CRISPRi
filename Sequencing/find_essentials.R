organism_hogs <- orthologs %>%
  select(HOG, locus_tag, organism) %>%
  inner_join(rbind(eco_results_median, ecl_results_median, kpn_results_median)) %>%
  data.table()

# try to find the E. coli "gene" (not locus_tag, select "gene") equivalent for each HOG to name the HOGs

e_coli_names <- organism_hogs %>%
  filter(organism == "E. coli") %>%
  select(gene, HOG) %>%
  unique() %>%
  rename(e_coli_name = gene) %>%
  # then collapse the gene names into a single comma separated list for each HOG
  group_by(HOG) %>%
  summarise(
    e_coli_name = list(sort(unique(e_coli_name)))
  ) %>%
  mutate(
    e_coli_name_group = vapply(
      e_coli_name,
      paste,
      collapse = ",",
      FUN.VALUE = character(1)
    )
  )


# create a single column with a comma-separated list of genes for each organism, gene names are found in organism_hogs$gene
complete_orthologs <- complete_orthologs %>%
  inner_join(organism_hogs) %>%
  group_by(HOG) %>%
  summarise(
    gene = list(sort(unique(gene)))
  ) %>%
  mutate(
    gene_group = vapply(
      gene,
      paste,
      collapse = ",",
      FUN.VALUE = character(1)
    )
  )






most_variable_HOGs <- organism_hogs %>%
  group_by(HOG, contrast) %>%
  summarise(
    logFC_var = var(logFC),
    logFC_mean = mean(logFC),
    gene_count = n()
  ) %>%
  arrange(desc(logFC_var))

# also list the variance of the logFC for each HOG in each contrast
HOG_variances <- organism_hogs %>%
  group_by(HOG, contrast) %>%
  summarise(
    logFC_var = var(logFC),
    logFC_mean = mean(logFC),
    gene_count = n()
  ) %>%
  arrange(desc(logFC_var))

HOG_variances %>%
  inner_join(complete_orthologs) %>%
  data.table() %>%
  dcast(HOG + gene_group ~ contrast, value.var = "logFC_var") %>%
  arrange(desc(induced_imipenem_1x))

########

ortho_enrichments %>%
  filter(term == "GO:0033592") %>%
  inner_join(organism_hogs) %>%
  inner_join(e_coli_names) %>%
  filter(contrast != "intercept") %>%
  filter(type %like% "perfect") %>%
  ggplot(aes(
    x = logFC,
    y = -log10(adj.P.Val)
  )) +
  geom_point() +
  facet_grid(contrast ~ organism) +
  # draw labels if adj.P.Val is less than 0.05
  geom_text_repel(
    # data = . %>% filter(adj.P.Val < 0.05),
    aes(label = e_coli_name_group),
    box.padding = 0.5,
    point.padding = 0.5
  )
