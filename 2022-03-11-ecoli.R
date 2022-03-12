well_locations <- 
  CJ(
    well_row = c(1:12), 
    well_col = toupper(letters[1:8]))[
      , .(well = paste0(well_col, well_row))]


orgs <- c(
  "5111",
  "5142")


drugs <- c(
  "none",
  "ampicillin",
  "mecillinam",
  "novobiocin",
  "urea",
  "vancomycin")


org_key <- data.table(
  organism = factor(orgs, levels = orgs))


drug_key <- data.table(
  drug = factor(drugs, levels = drugs))


layout <- CJ(
  org = org_key$organism, 
  drug = drug_key$drug, 
  induced = c(TRUE, FALSE), 
  rep = c(1:4))


layout <- cbind(well_locations, layout)


expanded_layout <- dcast(layout, induced + rep ~ org + drug, value.var = "well")


drug_exp <- fread(
  '2022-03-11-ecoli.tsv', 
  header = TRUE,
  na.strings = "NA")


setnames(drug_exp, t(drug_exp)[,1])


setnames(
  drug_exp,
  "Time [s]",
  "well")


drug_exp <-
  drug_exp[grep("[A-z][0-9]{1,2}", well)]


drug_exp <- 
  melt(
    drug_exp, 
    id.vars = "well", 
    variable.name = "time", 
    value.name = "OD600", 
    na.rm = TRUE)


drug_exp <- layout[drug_exp, on = .(well)]