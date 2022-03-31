library(pacman)

p_load(
  data.table, 
  RSQLite, 
  RColorBrewer,
  growthcurver,
  ggplot2)

chem_gen_db <- dbConnect(RSQLite::SQLite(), "chem_gen.db")

Experiments <- data.table(dbReadTable(chem_gen_db, "Experiments"))

Ryan_Strains <- data.table(dbReadTable(chem_gen_db, "Ryan_Strains"))

Well_Stats <- data.table(dbReadTable(chem_gen_db, "Well_Stats"))

dbDisconnect(chem_gen_db)

Well_Stats[, Row := gsub("[0-9]{1,2}", "", Well)]

Well_Stats[, Column := as.numeric(gsub("[A-Z]", "", Well))]

setorder(Well_Stats, Column, Row)

Well_Stats[, Summary := paste(Media, Rep, Induced, Organism, Chemical, Dose)]

