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
