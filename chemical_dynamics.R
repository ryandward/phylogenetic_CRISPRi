typas_top_cor <- 
  melt(
    data.table(
      typas_cor, 
      keep.rownames = "condition_x"), 
    id.vars = "condition_x", 
    variable.name = "condition_y",
    value.name = "correlation")

typas_top_cor[, c("condition_x", "dose_x") := tstrsplit(condition_x, "-")]
typas_top_cor[, c("condition_y", "dose_y") := tstrsplit(condition_y, "-")]

typas_top_cor <-
  typas_top_cor[
    condition_groups[
      , .(
        condition_x = condition, 
        dose_x = dose, 
        group_x = group)], 
    on = .(condition_x, dose_x)]

typas_top_cor <-
  typas_top_cor[
    condition_groups[
      , .(
        condition_y = condition, 
        dose_y = dose, 
        group_y = group)], 
    on = .(condition_y, dose_y)]

setorder(typas_top_cor, -correlation) 

intra_group_floor <- 
  typas_top_cor[
    correlation != 1 
    & group_x != 4 
    & group_y != 4
    & group_x != group_y, max(correlation)]
