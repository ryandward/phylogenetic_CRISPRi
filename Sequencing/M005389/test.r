# require(pacman)

# data %>%
#   inner_join(targets) %>%
#   # filter(sample %in% c("ECimiRDW1", "ECimiRDW4")) %>%
#   group_by(sample) %>%
#   mutate(cpm = cpm(count)[, 1]) %>%
#   group_by(locus_tag) %>%
#   filter(cpm ==
#     max(cpm)) %>%
#   select(spacer) %>%
#   unique() %>%
#   ungroup() %>%
#   select(spacer) %>%
#   inner_join(results) %>%
#   group_by(contrast, locus_tag) %>%
#   mutate(FDR = ifelse(FDR == 1, 0.99999, FDR)) %>%
#   inner_join(definitions) %>%
#   filter(contrast %in% c("induction_only", "imipenem_partial")) %>%
#   arrange(contrast, logFC) %>%
#   group_by(contrast) %>%
#   mutate(index_asc = row_number()) %>%
#   arrange(contrast, desc(logFC)) %>%
#   mutate(index_desc = row_number()) %>%
#   ungroup() %>%
#   mutate(gene_name = ifelse(is.na(gene_name), sprintf("bold('%s')", locus_tag), sprintf("italic('%s')", gene_name))) %>%
#   mutate(type = ifelse(type %in% c("mismatch", "perfect essential"), "essential", "nonessential")) %>%
#   ggplot(aes(y = logFC, x = FDR)) +
#   geom_hline(yintercept = 0, linetype = "solid", color = "grey50", size = 1.5) +
#   scale_x_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
#   geom_point(aes(fill = type, color = type), size = 2, alpha = 0.5, shape = 16) +
#   scale_color_manual(values = c("nonessential" = "#1F78B4", "essential" = "#E31A1C")) +
#   scale_fill_manual(values = c("nonessential" = "#1F78B4", "essential" = "#E31A1C")) +
#   theme_bw() +
#   facet_wrap(~contrast) +
#   theme(strip.text.y = element_text(angle = 0)) +
#   geom_point(
#     data = . %>% filter(FDR < 0.01 & abs(logFC) >= 1 & (index_asc <= 20 | index_desc <= 20)),
#     aes(color = type),
#     size = 4,
#     shape = 21,
#     alpha = 0.75,
#     stroke = 1
#   ) +
#   geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 0.5) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
#   geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", size = 0.5) +
#   geom_label_repel(
#     data = . %>% filter(FDR < 0.01 & abs(logFC) >= 1 & (index_asc <= 20 | index_desc <= 20)),
#     aes(label = gene_name),
#     size = 3,
#     segment.size = 0.5,
#     segment.color = "black",
#     segment.alpha = 0.5,
#     box.padding = 0.5,
#     point.padding = 0.5,
#     force = 10,
#     nudge_x = 0,
#     nudge_y = 0,
#     parse = TRUE
#   ) +
#   labs(x = "Confidence (FDR)", y = "Relative Fitness Score (log2FC)", color = "Type", fill = "Type")

# enrichments <- fread("/home/ryandward/Downloads/enrichment.all.tsv")

# # convert column "matching proteins in your input (IDs)" from csv to list
# enrichments <- enrichments %>%
#   mutate(
#     protein = strsplit(
#       x = `matching proteins in your input (IDs)`,
#       split = ", "
#     )
#   ) %>%
#   select(-`matching proteins in your input (IDs)`) %>%
#   unnest(protein)

# # replace protein with everything after the first period, i.e. delete text before first period
# enrichments <- enrichments %>%
#   mutate(protein = str_replace(protein, ".*\\.", "")) %>%
#   rename(locus_tag = protein)

# STRING <- enrichments %>%
#   filter(`#category` == "STRING clusters") %>%
#   select(`term description`, locus_tag)

# # only keep the first 20 characters of the term description
# STRING <- STRING %>%
#   mutate(`term description` = str_sub(`term description`, 1, 30))


# all_string <- fread("Organisms/511145.protein.enrichment.terms.v12.0.txt")

# enrichments <- all_string %>%
#   mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
#   filter(term %in% (all_string %>% group_by(term) %>% tally() %>% filter(n < 5) %>% pull(unique(term))))


# library(broom)
# library(purrr)

# aov_results <- data %>%
#   group_by(sample, target) %>%
#   mutate(count = sum(count)) %>%
#   ungroup() %>%
#   group_by(sample) %>%
#   mutate(cpm = cpm(count)) %>%
#   left_join(targets) %>%
#   right_join(enrichments) %>%
#   data.table() %>%
#   mutate(replicate = as.factor(replicate), imipenem = as.factor(imipenem), induced = as.factor(induced)) %>%
#   group_by(term) %>%
#   nest() %>%
#   mutate(aov_results = map(data, ~ {
#     aov_res <- aov(cpm ~ replicate + induced + imipenem, data = .)
#     results <- tidy(aov_res)
#     rename(results, aov_term = term)
#   })) %>%
#   select(-data) %>%
#   unnest(aov_results)

# good_terms <- aov_results %>%
#   filter(term %in% c(aov_results %>%
#     filter(aov_term == "replicate") %>%
#     filter(p.value > 0.10) %>%
#     arrange(p.value) %>% head(10) %>% unique() %>% pull(term))) %>%
#   pull(term) %>%
#   unique() %>%
#   c("PF00358")



# LRT_model <- data %>%
#   group_by(sample, target) %>%
#   mutate(count = sum(count)) %>%
#   ungroup() %>%
#   group_by(sample) %>%
#   mutate(cpm = cpm(count)) %>%
#   left_join(targets) %>%
#   right_join(enrichments) %>%
#   data.table() %>%
#   mutate(replicate = as.factor(replicate), imipenem = as.factor(imipenem), induced = as.factor(induced)) %>%
#   group_by(term) %>%
#   nest() %>%
#   mutate(aov_results = map(data, ~ {
#     reduced_model <- aov(cpm ~ replicate + induced + imipenem, data = .)
#     full_model <- aov(cpm ~ induced + imipenem, data = .)
#     lrt_result <- anova(reduced_model, full_model)
#     results <- tidy(lrt_result)
#     rename(results, aov_term = term)
#   })) %>%
#   select(-data) %>%
#   unnest(aov_results) %>%
#   filter(p.value <= 0.05)

# ###############

# library(purrr)

# all_string <- fread("Organisms/511145.protein.enrichment.terms.v12.0.txt") %>%
#   mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
#   unique()


# enrichments <- all_string %>%
#   filter(term %in% (all_string %>% group_by(term) %>% tally() %>% filter(n <= 100) %>% pull(unique(term)))) %>%
#   group_by(category, term, description) %>%
#   summarise(locus_tag = list(sort(unique(locus_tag)))) %>%
#   mutate(locus_tag_group = vapply(locus_tag, paste, collapse = ",", FUN.VALUE = character(1))) %>%
#   ungroup() %>%
#   distinct(locus_tag_group, .keep_all = TRUE) %>%
#   select(-locus_tag_group) %>%
#   unnest(locus_tag)

# enrichments_with_enough_spacers <- enrichments %>%
#   inner_join(targets, relationship = "many-to-many") %>%
#   unique() %>%
#   group_by(term) %>%
#   tally() %>%
#   filter(n >= 4) %>%
#   filter(!term %in% (enrichments %>% full_join(targets, relationship = "many-to-many") %>% filter(is.na(spacer)) %>% pull(term) %>% unique())) %>%
#   unique() %>%
#   arrange(n) %>%
#   inner_join(enrichments %>% select(term, description) %>% unique())

# safe_aov <- possibly(function(data) {
#   reduced_model <- aov(lcpm ~ induced, data = data)
#   full_model <- aov(lcpm ~ induced:imipenem, data = data)
#   lrt_result <- anova(reduced_model, full_model)
#   results <- broom::tidy(lrt_result)
#   rename(results, aov_term = term)
# }, otherwise = NULL)

# LRT_model <- full_data %>%
#   group_by(sample) %>%
#   mutate(lcpm = log(cpm(count+1), doublings)) %>%
#   inner_join(targets %>% select(spacer, locus_tag) %>% unique()) %>%
#   inner_join(enrichments) %>%
#   inner_join(enrichments_with_enough_spacers) %>%
#   data.table() %>%
#   mutate(replicate = as.factor(replicate), imipenem = as.factor(imipenem), induced = as.factor(induced)) %>%
#   group_by(term) %>%
#   nest() %>%
#   mutate(aov_results = map(data, safe_aov)) %>%
#   select(-data) %>%
#   unnest(aov_results) %>%
#   arrange(p.value) %>%
#   inner_join(enrichments %>% group_by(term, description) %>% tally(name = "gene_count"))

# good_terms <- NULL

# good_terms <- LRT_model %>%
#   filter(p.value <= 0.05 & aov_term == "lcpm ~ induced:imipenem") %>%
#   inner_join(enrichments_with_enough_spacers) %>%
#   group_by(description) %>%
#   slice_min(p.value, n = 1) %>%
#   ungroup() %>%
#   arrange(p.value) %>%
#   unique() %>%
#   head(35)


# good_terms$description <- factor(good_terms$description, levels = good_terms$description %>% unique())


# # Get the unique terms
# unique_terms <- unique(enrichments_with_enough_spacers$term)

# target_spacers_for_terms <- enrichments_with_enough_spacers %>%
#   inner_join(enrichments) %>%
#   inner_join(targets)

# # Initialize a list to store the gene indices
# gene_indices <- list()

# # Loop over the unique terms
# for (current_term in unique_terms) {
#   # Get the locus tags of the genes in the set to be tested
#   locus_tags <- target_spacers_for_terms %>%
#     filter(term == current_term) %>%
#     pull(spacer)

#   # Create a vector of indices for the genes in the set
#   gene_indices[[current_term]] <- which(rownames(dge) %in% locus_tags)
# }

# # Perform the competitive gene set test for all gene sets
# imipenem_partial_sets <- camera(dge, index = gene_indices, design = design_matrix, contrast = contrasts[, "imipenem_partial"]) %>%
#   data.table(keep.rownames = "term") %>%
#   mutate(term = factor(term, levels = unique_terms))

# induction_only_sets <- camera(dge, index = gene_indices, design = design_matrix, contrast = contrasts[, "induction_only"]) %>%
#   data.table(keep.rownames = "term") %>%
#   mutate(term = factor(term, levels = unique_terms))

# full_data %>%
#   group_by(sample) %>%
#   mutate(cpm = cpm(count + 1), doublings) %>%
#   ungroup() %>%
#   inner_join(targets) %>%
#   inner_join(enrichments) %>%
#   inner_join(imipenem_partial_sets %>% head(30) %>% select(term, FDR)) %>%
#   ggplot(aes(y = log(cpm, doublings), x = factor(induced))) +
#   geom_boxplot(aes(fill = factor(imipenem)), alpha = 0.5) +
#   facet_wrap(~ paste(term, paste0("(-logFDR=", round(-log10(FDR), 1), ")")) + stringr::str_sub(description, start = 1, end = 25), ncol = 6, scales = "free_y")


# full_data %>%
#   group_by(sample) %>%
#   mutate(cpm = cpm(count + 1), doublings) %>%
#   ungroup() %>%
#   inner_join(targets) %>%
#   inner_join(enrichments) %>%
#   inner_join(
#     induction_only_sets %>%
#       inner_join(imipenem_partial_sets, by = "term") %>%
#       mutate(diff = -log10(FDR.y) - -log10(FDR.x)) %>%
#       arrange(desc(diff)) %>%
#       inner_join(enrichments_with_enough_spacers) %>%
#       head(30) %>%
#       mutate(FDR = 10^(-diff))
#   ) %>%
#   arrange(FDR) %>%
#   mutate(facet_title = paste(term, paste0("(-logFDR=", round(-log10(FDR), 1), ")"))) %>%
#   mutate(facet_title = factor(facet_title, levels = facet_title %>% unique())) %>%
#   ggplot(aes(y = log(cpm, doublings), x = factor(induced))) +
#   geom_boxplot(aes(fill = factor(imipenem)), alpha = 0.5) +
#   facet_wrap(~ facet_title + stringr::str_sub(description, start = 1, end = 30), ncol = 6, scales = "free_y")


# results %>%
#   select(-PValue, -FDR) %>%
#   inner_join(targets) %>%
#   inner_join(enrichments) %>%
#   inner_join(imipenem_partial_sets %>% head(10) %>% select(term, FDR)) %>%
#   ggplot(aes(y = logFC)) +
#   geom_boxplot(aes(fill = contrast), alpha = 0.5) +
#   facet_wrap(~ term + description + paste("FDR =", FDR))




# Split the spacer column by term
locus_tag_lists <- split(targets$spacer, targets$locus_tag)

# Find the indices of each set of locus tags in rownames(dge)
locus_tag_indices <- lapply(locus_tag_lists, function(locus_tags) which(rownames(dge) %in% locus_tags))

locus_tag_v_fit <- voomWithQualityWeights(dge, design_matrix, plot=TRUE)

v_targets <- locus_tag_v_fit$E %>%
  data.table(keep.rownames = "spacer") %>%
  select(spacer) %>%
  left_join(
    targets %>% 
    filter(locus_tag %in% all_string$locus_tag) %>%
    group_by(spacer) %>%
    filter(target == "None" | (sp_dir!=tar_dir & abs(as.numeric(offset)) == min(abs(as.numeric(offset))) & overlap == max(overlap)))
    %>%
    group_by(target)) %>%
    ungroup %>%
    data.table()

v_targets[y_pred == "None", y_pred := NA_integer_]

v_targets$y_pred <- as.numeric(v_targets$y_pred)

v_targets[is.na(target) | target == "None", weight := min(v_targets$y_pred, na.rm = TRUE)]

v_targets[spacer == target, weight := max(v_targets$y_pred, na.rm = TRUE)]

v_targets[type == "mismatch", weight := y_pred]

v_targets$weight <- rescale(as.numeric(v_targets$weight), to = c(1, 100))



# Perform the competitive gene set test for all gene sets
locus_tag_confidence <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    locus_tag_v_fit, index = locus_tag_indices, design = design_matrix,
    # geneid = v_targets$locus_tag,
    # allow.neg.cor = TRUE,
    weights = v_targets$weight,
    # inter.gene.cor=0.05,
    # approx.zscore = FALSE,
    # trend.var = TRUE,
    # use.ranks = TRUE,
    contrast = contrast_column
  ) %>%
  data.table(keep.rownames = "locus_tag") %>% 
  mutate(contrast = factor(contrast_name))
  result
}) %>% 
do.call(rbind, .)

locus_tag_confidence <- v_targets %>% 
  select(locus_tag, weight, spacer) %>% unique() %>% 
  group_by(locus_tag) %>% 
  summarise(weights = sum(weight)) %>% 
  inner_join(locus_tag_confidence) %>% 
  rename(NGuides = NGenes) %>%
  data.table()


locus_tag_confidence %>% 
  arrange(FDR) %>% 
  inner_join(targets %>% select(locus_tag, gene) %>% unique()) %>% 
  inner_join(results %>% group_by(contrast, locus_tag) %>% 
  summarize(logFC = median(logFC))) %>% 
  mutate(gene = ifelse(is.na(gene), sprintf("bold('%s')", locus_tag), sprintf("italic('%s')", gene))) %>%
  mutate(type = ifelse(locus_tag %in% essentials, "essential", "nonessential")) %>%
  group_by(contrast) %>%
  arrange(contrast, logFC) %>%
  mutate(index_asc = row_number()) %>%
  arrange(contrast, desc(logFC)) %>%
  mutate(index_desc = row_number()) %>%
  ggplot(aes(y = logFC, x = FDR)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50", size = 1.5) +
  scale_x_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
  geom_point(aes(fill = type, color = type), size = 2, alpha = 0.5, shape = 16) +
  scale_color_manual(values = c("nonessential" = "#1F78B4", "essential" = "#E31A1C")) +
  scale_fill_manual(values = c("nonessential" = "#1F78B4", "essential" = "#E31A1C")) +
  theme_bw() +
  facet_wrap(~contrast, scales = "free") +
  theme(strip.text.y = element_text(angle = 0)) +
  geom_point(
    data = . %>% filter(FDR < 0.01 & abs(logFC) >= 1 & (index_asc <= 20 | index_desc <= 20)),
    aes(color = type),
    size = 4,
    shape = 21,
    alpha = 0.75,
    stroke = 1
  ) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", size = 0.5) +
  geom_label_repel(
    data = . %>% filter(FDR < 0.01 & abs(logFC) >= 1 & (index_asc <= 20 | index_desc <= 20)),
    aes(label = gene),
    size = 3,
    segment.size = 0.5,
    segment.color = "black",
    segment.alpha = 0.5,
    box.padding = 0.5,
    point.padding = 0.5,
    force = 10,
    nudge_x = 0,
    nudge_y = 0,
    parse = TRUE
  ) +
  labs(x = "Confidence (adj.P.Val)", y = "Relative Fitness Score (log2FC)", color = "Type", fill = "Type")





# enrichments %>%
#   inner_join(targets) %>%
#   filter(gene %like% "rpl") %>%
#   select(term) %>%
#   unique() %>%
#   inner_join(
#     all_sets %>%
#     filter(FDR <= 0.05) %>%
#     select(term, FDR, contrast)) %>%
#     inner_join(enrichments %>% select(-category)
#     ) %>%
#     inner_join(
#       results %>% select(logFC, locus_tag, spacer, adj.P.Val, contrast)
#       ) %>% 
#     group_by(
#       term, description, contrast, FDR
#       ) %>% 
#       summarise(medlogFC = median(logFC), sdlogFC = sd(logFC), n = n()) %>% arrange(sdlogFC) %>%
#       ggplot(aes(x = sdlogFC, y = -log10(FDR))) +
#       geom_point(
#         aes(size = n, fill = medlogFC),
#         color = "black", shape = 21) +
#       scale_size(range = c(2, 10)) +
#       facet_wrap(~contrast) + scale_fill_viridis() +
#       geom_label_repel(aes(label = description))


#########################

# edges <- enrichments %>%
#       inner_join(v_targets) %>%
#       # filter(gene %like% "proQ") %>%
#       select(term) %>%
#       unique() %>%
#       inner_join(
#         all_sets %>%
#         filter(FDR <= 0.05) %>%
#         select(term, FDR, contrast)) %>%
#         inner_join(enrichments %>% select(-category)
#         ) %>%
#         inner_join(
#           results %>% select(logFC, locus_tag, gene, spacer, adj.P.Val, contrast) %>%
#           inner_join(v_targets %>% select(spacer, weight) %>% unique())
#           ) %>% 
#           filter(adj.P.Val <= 0.05) %>%
#           filter(contrast == "induced")


edges <- (all_sets[contrast == "induced"] %>% 
  rename(FDR.term = FDR, PValue.term = PValue, Direction.term = Direction)) %>% 
  inner_join(enrichments %>% select(term, description, category, locus_tag) %>% unique()) %>%
  inner_join(locus_tag_confidence[contrast == "induced"] %>% 
  rename(FDR.gene= FDR, PValue.gene = PValue, Direction.gene = Direction)) %>% 
  filter(FDR.term <= 0.05 & FDR.gene <= 0.05) %>%
  inner_join(v_targets %>% select(locus_tag, gene) %>% 
  unique()) %>%
  filter(Direction.term == Direction.gene) 
  
terms <- unique(edges$term)
genes <- unique(edges$gene)

# Create an edge list from your data
edge_list <- edges %>% select(term, gene, weights)

# Create the graph
g <- graph_from_data_frame(edge_list, directed = FALSE)

# Assign node types based on whether they are in the terms or locus_tags list
V(g)$type <- ifelse(V(g)$name %in% terms, "term", "gene")

# Assign attributes to nodes
V(g)$FDR <- ifelse(V(g)$type == "term", edges$FDR.term[match(V(g)$name, edges$term)], edges$FDR.gene[match(V(g)$name, edges$gene)])
# V(g)$logFC <- ifelse(V(g)$type == "term", NA, edges$logFC[match(V(g)$name, edges$gene)])

# Plot the graph
plot(g, vertex.size = 2, layout = layout_with_kk, vertex.color = ifelse(V(g)$type == "term", "lightgreen", "lightblue"))


# Identify the connected components
comp <- components(g)

# Split the graph into a list of subgraphs and calculate centrality measures
subgraphs <- lapply(1:max(comp$membership), function(i) {
  subgraph <- induced_subgraph(g, which(comp$membership == i))
  # Calculate degree centrality
  V(subgraph)$degree <- degree(subgraph)
  # Calculate closeness centrality
  V(subgraph)$closeness <- closeness(subgraph)
  # Calculate betweenness centrality
  V(subgraph)$betweenness <- betweenness(subgraph)
  # Calculate eigenvector centrality
  V(subgraph)$eigen_centrality <- eigen_centrality(subgraph)$vector
  return(subgraph)
})

# Extract all vertex attributes for each subgraph and combine them into a long-form data.table
vertex_attr_dt <- rbindlist(lapply(1:length(subgraphs), function(i) {
  dt <- as.data.table(vertex_attr(subgraphs[[i]]))
  dt[, graph := i]
  return(dt)
}), fill = TRUE)

# vertex_attr_dt %>% arrange(desc(betweenness)) %>% filter(type == "term") %>% rename(term = name) %>% inner_join(all_string %>% group_by(description, term, c
#     ategory) %>% tally()) %>% filter(degree == n)



#########################

# Create a term-term matrix
gene_gene_matrix <- table(edges$gene, edges$term)
term_term_matrix <- table(edges$term, edges$gene)

# Convert the term-term matrix to a term-term adjacency matrix
gene_adj_matrix <- gene_gene_matrix %*% t(gene_gene_matrix)
term_adj_matrix <- term_term_matrix %*% t(term_term_matrix)

# Create a graph from the adjacency matrix
gene_adj_g <- graph_from_adjacency_matrix(gene_adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
term_adj_g <- graph_from_adjacency_matrix(term_adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Plot the graph
plot(gene_adj_g, vertex.size = 2, layout = layout_with_kk, vertex.color = "lightgreen")
plot(term_adj_g, vertex.size = 2, layout = layout_with_kk, vertex.color = "lightgreen")

# Identify the connected components
gene_adj_comp <- components(gene_adj_g)
term_adj_comp <- components(term_adj_g)


gene_adj_subgraphs <- lapply(1:max(gene_adj_comp$membership), function(i) {
  gene_adj_subgraph <- induced_subgraph(gene_adj_g, which(gene_adj_comp$membership == i))
  # Calculate degree centrality
  V(gene_adj_subgraph)$degree <- degree(gene_adj_subgraph)
  # Calculate closeness centrality
  V(gene_adj_subgraph)$closeness <- closeness(gene_adj_subgraph)
  # Calculate betweenness centrality
  V(gene_adj_subgraph)$betweenness <- betweenness(gene_adj_subgraph)
  # Calculate eigenvector centrality
  V(gene_adj_subgraph)$eigen_centrality <- eigen_centrality(gene_adj_subgraph)$vector
  return(gene_adj_subgraph)
})

# Split the graph into a list of subgraphs and calculate centrality measures
term_adj_subgraphs <- lapply(1:max(term_adj_comp$membership), function(i) {

  term_adj_subgraph <- induced_subgraph(term_adj_g, which(term_adj_comp$membership == i))
  # Calculate degree centrality
  V(term_adj_subgraph)$degree <- degree(term_adj_subgraph)
  # Calculate closeness centrality
  V(term_adj_subgraph)$closeness <- closeness(term_adj_subgraph)
  # Calculate betweenness centrality
  V(term_adj_subgraph)$betweenness <- betweenness(term_adj_subgraph)
  # Calculate eigenvector centrality
  V(term_adj_subgraph)$eigen_centrality <- eigen_centrality(term_adj_subgraph)$vector
  return(term_adj_subgraph)
})

gene_adj_verts <- rbindlist(lapply(seq_along(gene_adj_subgraphs), function(i) {
  dt <- as.data.table(vertex_attr(gene_adj_subgraphs[[i]]))
  dt[, graph := i]
  return(dt)
}), fill = TRUE)

term_adj_verts <- rbindlist(lapply(seq_along(term_adj_subgraphs), function(i) {
  dt <- as.data.table(vertex_attr(term_adj_subgraphs[[i]]))
  dt[, graph := i]
  return(dt)
}), fill = TRUE)

 
term_adj_summary <- vertex_attr(term_adj_g) %>% 
  list() %>% 
  rbindlist() %>% 
  arrange(desc(closeness)) %>% 
  rename(term = name) %>% 
  inner_join(enrichments %>% 
               group_by(term, description) %>% 
               tally()) %>% 
  inner_join(enrichments %>% 
               inner_join(v_targets %>% 
                            select(locus_tag, gene) %>% 
                            unique())) %>%
  data.table() %>%
  .[, locus_tags := paste(locus_tag, collapse = ","), by = term] %>%
  .[, genes := paste(gene, collapse = ","), by = term] %>%
  .[, c("locus_tag", "gene") := NULL] %>%
  unique() %>%
  inner_join(all_sets[contrast == "induced"])





