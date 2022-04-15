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

Experiments.pk <- c("Date", "Well", "Time")

Experiments[, Organism := as.integer(Organism)]

Well_Stats <- 
  Experiments[
    , SummarizeGrowth(
      Time/60/60, 
      OD600)$vals, 
    by = .(
      Date,
      Well,
      cJMP,
      Chemical,
      Dose,
      Unit,
      Instrument,
      Rep,
      Media,
      Induced,
      Organism)]


Fitted_Experiments <- 
  Well_Stats[
    , .(
      Hour = c(1:15), 
      OD600_fit = k  / ( 1 + ( ( k - n0 ) / n0 ) * exp( -r * c(1:15) ) ) ), 
    by= .(
      Date,
      Well,
      cJMP,
      Chemical,
      Dose,
      Unit,
      Instrument,
      Rep,
      Media,
      Induced,
      Organism)]

Fitted_Experiments[
  is.na(OD600_fit), 
  OD600_fit := 0]

dbWriteTable(
  chem_gen_db,
  "Well_Stats",
  Well_Stats,
  overwrite = TRUE)

dbWriteTable(
  chem_gen_db,
  "Fitted_Experiments",
  Fitted_Experiments,
  overwrite = TRUE)

dbDisconnect(chem_gen_db)