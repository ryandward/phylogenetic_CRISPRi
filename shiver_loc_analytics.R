shiver <- protein_loc[shiver, on = .(gene == Gene)]
setnames(shiver, "gene", "Gene")

significant_conditions_loc <-
  shiver[abs(Score) > (2 * sd(shiver$Score, na.rm = T)), 
             .N, 
             by = .(Condition, STEPdb_loc_code)]

significant_conditions_loc <- 
  significant_conditions_loc[condition_groups, on = .(Condition)]

most_phenotypes_by_group_loc <- 
  significant_conditions_loc[significant_conditions_loc[, .(N = max(N, na.rm = TRUE)), by = .(Group, STEPdb_loc_code)], on = .(N, STEPdb_loc_code, Group)]

#https://docs.google.com/document/d/1hvRGJ7ZJy1We3j9Vh1qcUzapQK17SXealcnLI-lOzBw/edit#