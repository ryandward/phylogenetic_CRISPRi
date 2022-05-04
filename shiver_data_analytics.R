library('pacman')

p_load(
  data.table, 
  pheatmap,
  factoextra)


#https://pubmed.ncbi.nlm.nih.gov/27355376/
shiver <- fread("shiver.tsv.gz")

################################################################################
# how many groups for k-means? 
target_k = 10
################################################################################


shiver <- melt(
  shiver, 
  variable.name = "Gene", 
  value.name = "Score")


shiver[, c("Gene", "Other") := tstrsplit(Gene, "\\{")]
shiver[, c("Gene", "Other") := tstrsplit(Gene, "-")]
shiver[, c("Gene", "Other") := tstrsplit(Gene, "\\*")]


shiver_df <- dcast(
  shiver, 
  Gene ~ Condition, 
  value.var = "Score", 
  fun.aggregate = median)


shiver_df_matrix <- data.matrix(shiver_df[, -1])


rownames(shiver_df_matrix) <- shiver_df$Gene


# shiver_df_matrix <- t(t(shiver_df_matrix)[complete.cases(t(shiver_df_matrix)),])


shiver_cor <- cor(shiver_df_matrix, use = "pairwise.complete.obs")


plot_matrix <- shiver_cor


to_plot <- pheatmap(
  plot_matrix,
  clustering_method = "ward.D2",
  clustering_distance_rows = "maximum",
  clustering_distance_cols = "maximum"
)


target_h <-
  sort(
    to_plot$tree_row$height, 
    decreasing = T)[target_k]

plot(
  to_plot$tree_row, 
  cex = 0.35,
  main = "Condition Similarity from Shiver Dataset")

abline(
  h = target_h, 
  col = "red", 
  lty = 2, 
  lwd = 2)

condition_groups <- 
  data.matrix(
    sort(
      cutree(
        to_plot$tree_row, 
        h = target_h)))

condition_groups <- 
  data.table(
    condition_groups, 
    keep.rownames = "Condition")

setnames(
  condition_groups, 
  "V1", 
  "Group")



significant_conditions <- 
  shiver[
    abs(Score) > (1 * sd(shiver$Score, na.rm = T)), 
    .N, 
    by = .(Condition)]


setorder(significant_conditions, N)


significant_conditions <- 
  significant_conditions[condition_groups, on = .(Condition)]

significant_conditions <- 
  significant_conditions[, .(Condition, N, Group)]

most_phenotypes_by_group <-
  significant_conditions[significant_conditions[, .(N = max(N, na.rm = TRUE)), by = .(Group)], on = .(Group, N)]




