today = "2022-03-09"
organism = "ecoli"

drug_exp <- fread(
  paste0(today, "-", organism, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

orgs <- c(
  "5111",
  "5142")

drugs <- c(
  "none",
  "deoxycholate",
  "mecillinam",
  "nickel(ii)",
  "novobiocin",
  "urea")

concs <- c(
  NA,
  1, 
  30, 
  1, 
  30, 
  750)

units <- c(
  NA,
  "%w/v",
  "ng/mL",
  "mM",
  "ug/mL",
  "mM")

org_key <- data.table(
  organism = factor(orgs, levels = unique(orgs)))

drug_key <- data.table(
  drug = factor(drugs, levels = unique(drugs)))

conc_key <- data.table(
  conc = factor(concs, levels = unique(concs)))

unit_key <- data.table(
  unit = factor(units, levels = unique(units)))

drug_key <- cbind(drug_key, conc_key, unit_key)

well_locations <- 
  CJ(
    well_row = c(1:12), 
    well_col = toupper(letters[1:8]))[
      , .(well = paste0(well_col, well_row))]

layout <- CJ(
  org = org_key$organism, 
  drug = drug_key[, paste(drug, conc, unit, sep = "_")], 
  induced = c(TRUE, FALSE), 
  rep = c(1:4))

layout <- cbind(well_locations, layout)

expanded_layout <- dcast(layout, induced + rep ~ org + drug, value.var = "well")

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

drug_exp[, c("drug", "dose", "units") := tstrsplit(drug, "_")]

drug_exp[, date := today]

drug_exp <- drug_exp[, .(date, well, org, drug, dose, units, induced, rep, time, OD600)]