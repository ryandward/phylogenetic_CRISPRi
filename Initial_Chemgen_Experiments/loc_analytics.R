################################################################################
# how many groups for k-means? 
target_k = 10
################################################################################

protein_dict <- fread("protein_dict.tsv")
protein_loc <- fread("protein_localization.tsv")
protein_loc <- protein_loc[protein_dict, on = .(protein)]

typas_OMP_OMLP <- typas_data[STEPdb_loc_code %in% c("OMLP", "OMP")]

typas_df_OMP_OMLP <- dcast(typas_OMP_OMLP, b_name ~ condition + dose, value.var = "score")

typas_df_OMP_OMLP_matrix <- data.matrix(typas_df_OMP_OMLP[, -1])

rownames(typas_df_OMP_OMLP_matrix) <- typas_df_OMP_OMLP$b_name

typas_df_OMP_OMLP_matrix <- typas_df_OMP_OMLP_matrix[complete.cases(typas_df_OMP_OMLP_matrix),]

typas_cor <- cor(typas_df_OMP_OMLP_matrix, use = "complete.obs")

plot_matrix <- typas_cor

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
  cex = 0.35)

abline(
  h = target_h, 
  col="red", 
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
    keep.rownames = "condition")

setnames(
  condition_groups, 
  "V1", 
  "group")

condition_groups[, c("condition", "dose") := tstrsplit(condition, "_")]

# fwrite(
#   condition_groups, 
#   "condition_groups.tsv", 
#   sep = "\t")
# 


# typas_data <- protein_loc[typas_data, on = .(gene == long_name)]

significant_conditions_OMP_OMLP <- 
  typas_OMP_OMLP[
    abs(score) > (1 * sd(typas_data$score, na.rm = T)), 
    .N, 
    by = .(condition, dose)]


setorder(significant_conditions_OMP_OMLP, N)


significant_conditions_OMP_OMLP <- 
  significant_conditions_OMP_OMLP[condition_groups, on = .(condition, dose)]

significant_conditions_OMP_OMLP <- 
  significant_conditions_OMP_OMLP[, .(condition, dose, N, group)]

most_phenotypes_by_group_OMP_OMLP <-
  significant_conditions_OMP_OMLP[significant_conditions_OMP_OMLP[, .(N = max(N, na.rm = TRUE)), by = .(group)], on = .(group, N)]

