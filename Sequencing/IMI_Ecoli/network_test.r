library(igraph)

# Split the spacer column by term
locus_tag_lists <- split(targets$spacer, targets$locus_tag)

# Find the indices of each set of locus tags in rownames(dge)
locus_tag_indices <- lapply(locus_tag_lists, function(locus_tags) which(rownames(dge) %in% locus_tags))

# locus_tag_v_fit <- voomWithQualityWeights(dge, design_matrix, plot=TRUE)

# v_targets <- locus_tag_v_fit$E %>%
#   data.table(keep.rownames = "spacer") %>%
#   select(spacer) %>%
#   left_join(
#     targets %>% 
#     filter(locus_tag %in% all_string$locus_tag) %>%
#     group_by(spacer) %>%
#     filter(target == "None" | (sp_dir!=tar_dir & abs(as.numeric(offset)) == min(abs(as.numeric(offset))) & overlap == max(overlap)))
#     %>%
#     group_by(target)) %>%
#     ungroup %>%
#     data.table()

# v_targets[y_pred == "None", y_pred := NA_integer_]

# v_targets$y_pred <- as.numeric(v_targets$y_pred)

# v_targets[is.na(target) | target == "None", weight := min(v_targets$y_pred, na.rm = TRUE)]

# v_targets[spacer == target, weight := max(v_targets$y_pred, na.rm = TRUE)]

# v_targets[type == "mismatch", weight := y_pred]

# v_targets$weight <- rescale(as.numeric(v_targets$weight), to = c(1, 100))



# Perform the competitive gene set test for all gene sets
locus_tag_confidence <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v, index = locus_tag_indices, design = design_matrix,
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
  rename(NGuides_gene = NGenes) %>%
  data.table()


# locus_tag_confidence %>% 
#   arrange(FDR) %>% 
#   inner_join(targets %>% select(locus_tag, gene) %>% unique()) %>% 
#   inner_join(results %>% group_by(contrast, locus_tag) %>% 
#   summarize(logFC = median(logFC))) %>% 
#   mutate(gene = ifelse(is.na(gene), sprintf("bold('%s')", locus_tag), sprintf("italic('%s')", gene))) %>%
#   mutate(type = ifelse(locus_tag %in% essentials, "essential", "nonessential")) %>%
#   group_by(contrast) %>%
#   arrange(contrast, logFC) %>%
#   mutate(index_asc = row_number()) %>%
#   arrange(contrast, desc(logFC)) %>%
#   mutate(index_desc = row_number()) %>%
#   ggplot(aes(y = logFC, x = FDR)) +
#   geom_hline(yintercept = 0, linetype = "solid", color = "grey50", size = 1.5) +
#   scale_x_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
#   geom_point(aes(fill = type, color = type), size = 2, alpha = 0.5, shape = 16) +
#   scale_color_manual(values = c("nonessential" = "#1F78B4", "essential" = "#E31A1C")) +
#   scale_fill_manual(values = c("nonessential" = "#1F78B4", "essential" = "#E31A1C")) +
#   theme_bw() +
#   facet_wrap(~contrast, scales = "free") +
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
#     aes(label = gene),
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
#   labs(x = "Confidence (adj.P.Val)", y = "Relative Fitness Score (log2FC)", color = "Type", fill = "Type")





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
#           filter(contrast == "induced_imipenem")


edges <- (all_sets[contrast == "induced_imipenem"] %>% 
  rename(FDR.term = FDR, PValue.term = PValue, Direction.term = Direction)) %>% 
  inner_join(enrichments %>% select(term, description, category, locus_tag) %>% unique()) %>%
  inner_join(locus_tag_confidence[contrast == "induced_imipenem"] %>% 
  rename(FDR.gene= FDR, PValue.gene = PValue, Direction.gene = Direction)) %>% 
  filter(FDR.term <= 0.01 & FDR.gene <= 0.01) %>%
  inner_join(v_targets %>% select(locus_tag, gene) %>% 
  unique()) %>%
  # filter(Direction.term == Direction.gene) %>%
  mutate(description = paste(description, term))

terms <- unique(edges$term)
genes <- unique(edges$gene)

# Create an edge list from your data
edge_list <- edges %>% filter (term %like% ("CL:")) %>% select(description, gene)
# edge_list %>% fwrite("edges.tsv", sep = "\t")

# Create the graph
g <- graph_from_data_frame(edge_list, directed = FALSE)

# Assign node types based on whether they are in the terms or locus_tags list
V(g)$type <- ifelse(V(g)$name %in% terms, "term", "gene")

# Assign attributes to nodes
V(g)$FDR <- ifelse(V(g)$type == "term", edges$FDR.term[match(V(g)$name, edges$term)], edges$FDR.gene[match(V(g)$name, edges$gene)])
# V(g)$logFC <- ifelse(V(g)$type == "term", NA, edges$logFC[match(V(g)$name, edges$gene)])

V(g)$FDR <- -log10(V(g)$FDR)

V(g)$Direction <- ifelse(V(g)$type == "term", edges$Direction.term[match(V(g)$name, edges$term)], edges$Direction.gene[match(V(g)$name, edges$gene)])

edge_list %>% fwrite("induced_imipenem_edges_CL.tsv", sep = "\t")

gene_nodes <- edges %>% 
  rename(name = gene, annotation = locus_tag, FDR = FDR.gene, Direction = Direction.gene) %>%
  select(name, annotation, FDR, Direction)

gene_nodes[, type := "gene"]

term_nodes <- edges %>%
  rename(name = description, annotation = term, FDR = FDR.term, Direction = Direction.term) %>%
  select(name, annotation, FDR, Direction)

term_nodes[, type := "term"]

gene_nodes %>% rbind(term_nodes) %>% mutate(FDR = -log10(FDR)) %>%
fwrite("induced_imipenem_node_info_CL.tsv", sep = "\t")

###################


write.graph(g, file = "temp_network.graphml", format = "graphml")

# Plot the graph
plot(g, 
vertex.size = 2, 
vertex.label = NA,
layout = layout_with_kk, 
vertex.color = ifelse(V(g)$type == "term", "lightgreen", "lightblue"))


# Identify the connected components
comp <- components(g)

# Split the graph into a list of subgraphs and calculate centrality measures
subgraphs <- lapply(1:max(comp$membership), function(i) {
  subgraph <- induced_imipenem_subgraph(g, which(comp$membership == i))
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
  gene_adj_subgraph <- induced_imipenem_subgraph(gene_adj_g, which(gene_adj_comp$membership == i))
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

  term_adj_subgraph <- induced_imipenem_subgraph(term_adj_g, which(term_adj_comp$membership == i))
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
  inner_join(all_sets[contrast == "induced_imipenem"])





