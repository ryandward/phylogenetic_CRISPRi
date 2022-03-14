library(pacman)
p_load(data.table, RSQLite)

chem_gen_db <- dbConnect(RSQLite::SQLite(), "chem_gen.db")

Experiments <- data.table(dbReadTable(chem_gen_db, "Experiments"))
Ryans_Chemicals <- data.table(dbReadTable(chem_gen_db, "Ryans_Chemicals"))


dbDisconnect(chem_gen_db)

Solution_Setup <- fread ('target_doses.tsv')
setorder(Solution_Setup, Chemical, Dose)

Organisms <- c("5147", "5142")

Chemicals <- Solution_Setup$Chemical
Chemicals <- c("none", Chemicals)

Doses <- Solution_Setup$Dose
Doses <- c("0", Doses)


Solution_Setup <- Ryans_Chemicals[, .(Chemical, Unit, cJMP, Stock_Concentration, Solvent)][Solution_Setup, on = .(Chemical)]

Solution_Setup[, Dilution := Stock_Concentration/Dose]

# Solution_Setup[, Max_Dilution := max(Dilution, na.rm = T), by = .(Date)]

Solution_Setup[, Work_Conc := 4]

Solution_Setup[, Work_Vol := 5000]

Solution_Setup[, Stock_Vol := Work_Conc * Work_Vol / Dilution]

Solution_Setup[, Extra_Water := Work_Vol - Stock_Vol]

Organism_key <- data.table(
  Organism = factor(Organisms, levels = unique(Organisms)))

Chemical_key <- data.table(
  Chemical = factor(Chemicals, levels = unique(Chemicals)))

Dose_key <- data.table(
  Dose = factor(Doses, levels = unique(Doses)))

Chemical_key <- cbind(Chemical_key, Dose_key)

Well_locations <- 
  CJ(
    Well_row = c(1:12), 
    Well_col = toupper(letters[1:8]))[
      , .(Well = paste0(Well_col, Well_row))]

layout <- CJ(
  Organism = Organism_key$Organism, 
  Chemical = factor(
    Chemical_key[, paste(Chemical, Dose, sep = ":")],
    levels = Chemical_key[, paste(Chemical, Dose, sep = ":")])
  , 
  Induced = c(TRUE, FALSE), 
  Rep = c(1:4))

layout <- cbind(Well_locations, layout)

expanded_layout <- dcast(layout, Induced + Rep ~ Organism + Chemical, value.var = "Well")
