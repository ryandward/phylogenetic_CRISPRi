typas_data <- protein_loc[typas_data, on = .(gene == long_name)]

significant_conditions_loc <-
  typas_data[abs(score) > (2 * sd(typas_data$score, na.rm = T)), 
             .N, 
             by = .(condition, STEPdb_loc_code)]

significant_conditions_loc <- 
  significant_conditions_loc[condition_groups, on = .(condition)]

most_phenotypes_by_group_loc <- 
  significant_conditions_loc[significant_conditions_loc[, .(N = max(N, na.rm = TRUE)), by = .(group, STEPdb_loc_code)], on = .(N, STEPdb_loc_code, group)]

