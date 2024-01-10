# Split the spacer column by term
locus_tag_lists <- split(targets$spacer, targets$locus_tag)

# Find the indices of each set of locus tags in rownames(dge)
locus_tag_indices <- lapply(locus_tag_lists, function(locus_tags) which(rownames(dge) %in% locus_tags))

# Perform the competitive gene set test for all gene sets
locus_tag_confidence <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v,
    index = locus_tag_indices, design = design_matrix,
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
  select(locus_tag, weight, spacer) %>%
  unique() %>%
  group_by(locus_tag) %>%
  summarise(weights = sum(weight)) %>%
  inner_join(locus_tag_confidence) %>%
  rename(NGuides = NGenes) %>%
  data.table()


# Get the unique contrasts
# contrasts <- unique(all_sets$contrast)

# Function to create a graph for a given contrast
create_graph <- function(current_contrast) {
  edges <- (all_sets[current_contrast == contrast] %>%
    rename(FDR.term = FDR, PValue.term = PValue, Direction.term = Direction)) %>%
    inner_join(enrichments %>% select(term, description, category, locus_tag) %>% unique()) %>%
    inner_join(locus_tag_confidence[current_contrast == contrast] %>%
      rename(FDR.gene = FDR, PValue.gene = PValue, Direction.gene = Direction)) %>%
    filter(FDR.term <= 0.05 & FDR.gene <= 0.05) %>%
    inner_join(v_targets %>% select(locus_tag, gene) %>%
      unique()) %>%
    filter(Direction.term == Direction.gene)

  # Create an edge list from your data
  edge_list <- edges %>% select(term, locus_tag, weights)

  # Create the graph
  g <- graph_from_data_frame(edge_list, directed = FALSE)

  # Add contrast as a graph attribute
  g$contrast <- current_contrast

  return(g)
}

# Create a list of graphs, one for each contrast
graph_list <- lapply(colnames(contrasts), create_graph)


# Function to assign node types and attributes
assign_attributes <- function(g) {
  # Assign node types
  V(g)$type <- ifelse(V(g)$name %in% enrichments$term, "term", "locus_tag")

  # Assign attributes to nodes
  # V(g)$FDR <- ifelse(V(g)$type == "term", edges$FDR.term[match(V(g)$name, edges$term)], edges$FDR.gene[match(V(g)$name, edges$gene)])
  # V(g)$logFC <- ifelse(V(g)$type == "term", NA, edges$logFC[match(V(g)$name, edges$gene)])

  return(g)
}

# Apply the function to each graph in the list
graph_list <- lapply(graph_list, assign_attributes)


# Function to calculate centrality measures and extract vertex attributes
calculate_centrality_and_extract_attributes <- function(g) {
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
    dt[, contrast := g$contrast] # Add the contrast column
    return(dt)
  }), fill = TRUE)

  return(list(subgraphs = subgraphs, vertex_attr_dt = vertex_attr_dt))
}

# Apply the function to each graph in the list
result_list <- lapply(graph_list, calculate_centrality_and_extract_attributes)

# Extract the vertex_attr_dt from each result in the list
vertex_attr_dt_list <- lapply(result_list, function(x) x$vertex_attr_dt)

# Combine all the data tables into one large data table
combined_vertex_attr_dt <- rbindlist(vertex_attr_dt_list, fill = TRUE)



create_adjacency_graphs <- function(g) {
  # Extract the edges from the graph
  edges <- as.data.frame(get.edgelist(g, names = TRUE))

  # Determine whether V1 or V2 is a gene or a term
  edges$V1_type <- V(g)$type[match(edges$V1, V(g)$name)]
  edges$V2_type <- V(g)$type[match(edges$V2, V(g)$name)]

  # Ensure that genes are in the first column and terms are in the second column
  edges_ordered <- edges
  if (unique(edges$V1_type) == "term") {
    edges_ordered <- edges[, c("V2", "V1", "V2_type", "V1_type")]
    names(edges_ordered) <- c("V1", "V2", "V1_type", "V2_type")
  }

  # Create a gene-term matrix
  gene_term_matrix <- table(edges_ordered$V1, edges_ordered$V2)

  # Convert the gene-term matrix to a gene-gene adjacency matrix and a term-term adjacency matrix
  gene_adj_matrix <- gene_term_matrix %*% t(gene_term_matrix)
  term_adj_matrix <- t(gene_term_matrix) %*% gene_term_matrix

  # Create a graph from the adjacency matrix
  gene_adj_g <- graph_from_adjacency_matrix(gene_adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  term_adj_g <- graph_from_adjacency_matrix(term_adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

  # Add the contrast attribute to the graphs
  graph_attr(gene_adj_g, "contrast") <- graph_attr(g, "contrast")
  graph_attr(term_adj_g, "contrast") <- graph_attr(g, "contrast")

  return(list(gene_adj_g = gene_adj_g, term_adj_g = term_adj_g))
}

# Apply the function to each graph in the list
adjacency_graph_list <- lapply(graph_list, create_adjacency_graphs)


perform_operations <- function(adjacency_graphs) {
  # Extract the adjacency graphs
  gene_adj_g <- adjacency_graphs$gene_adj_g
  term_adj_g <- adjacency_graphs$term_adj_g

  # Extract the contrast attribute
  contrast <- graph_attr(gene_adj_g, "contrast")

  # Identify the connected components
  gene_adj_comp <- components(gene_adj_g)
  term_adj_comp <- components(term_adj_g)

  # Split the graph into a list of subgraphs and calculate centrality measures
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

  # Combine all vertex attributes into a data.table
  gene_adj_verts <- rbindlist(lapply(seq_along(gene_adj_subgraphs), function(i) {
    dt <- as.data.table(vertex_attr(gene_adj_subgraphs[[i]]))
    dt[, graph := i]
    dt[, contrast := contrast]
    return(dt)
  }), fill = TRUE)

  term_adj_verts <- rbindlist(lapply(seq_along(term_adj_subgraphs), function(i) {
    dt <- as.data.table(vertex_attr(term_adj_subgraphs[[i]]))
    dt[, graph := i]
    dt[, contrast := contrast]
    return(dt)
  }), fill = TRUE)

  return(list(gene_adj_verts = gene_adj_verts, term_adj_verts = term_adj_verts))
}

# Apply the function to each pair of adjacency graphs in the list
result_list <- lapply(adjacency_graph_list, perform_operations)

##############################
# Extract all term_adj_verts and gene_adj_verts from the result list
term_adj_verts_list <- lapply(result_list, function(x) x$term_adj_verts)
gene_adj_verts_list <- lapply(result_list, function(x) x$gene_adj_verts)

# Combine all term_adj_verts and gene_adj_verts into a single data table
combined_term_adj_verts <- rbindlist(term_adj_verts_list, fill = TRUE)
combined_gene_adj_verts <- rbindlist(gene_adj_verts_list, fill = TRUE)


gene_centrality <- combined_gene_adj_verts %>%
  rename(locus_tag = name) %>%
  inner_join(targets %>% select(locus_tag, gene) %>% unique()) %>%
  select(contrast, locus_tag, gene, graph, degree, closeness, betweenness, eigen_centrality)

term_centrality <- combined_term_adj_verts %>%
  rename(term = name) %>%
  inner_join(enrichments %>% select(term, description) %>% unique()) %>%
  # mutate(missing_genes = gene_count - genes_targeted) %>%
  # arrange(missing_genes, desc(genes_targeted), desc(betweenness)) %>%
  inner_join(all_sets %>% select(term, contrast, FDR, Direction))

term_centrality %>%
  select(term, description, Direction, FDR, contrast, eigen_centrality, betweenness, graph) %>%
  as_tibble() %>%
  print(n = 25)
