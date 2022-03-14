library(pacman)
p_load(data.table, RSQLite)

chem_gen_db <- dbConnect(RSQLite::SQLite(), "chem_gen.db")

Experiments <- data.table(dbReadTable(chem_gen_db, "Experiments"))
Experiments.pk <- c("Date", "Well", "Time")
 
Experiments[, Organism := factor(Organism)]
# Experiments[, Induced := as.logical(Induced)]

Ryans_Chemicals <- data.table(dbReadTable(chem_gen_db, "Ryans_Chemicals"))

Organism_key <- data.table(
  Organism = factor(Organisms, levels = unique(Organisms)))

Chemical_key <- data.table(
  Chemical = factor(Chemicals, levels = unique(Chemicals)))

Dose_key <- data.table(
  Dose = factor(Doses, levels = unique(Doses)))

Unit_key <- data.table(
  Unit = factor(Units, levels = unique(Units)))

Chemical_key <- cbind(Chemical_key, Dose_key, Unit_key)

Well_locations <- 
  CJ(
    Well_row = c(1:12), 
    Well_col = toupper(letters[1:8]))[
      , .(Well = paste0(Well_col, Well_row))]

layout <- CJ(
  Organism = Organism_key$Organism, 
  Chemical = factor(
    Chemical_key[, paste(Chemical, Dose, Unit, sep = "_")],
    levels = Chemical_key[, paste(Chemical, Dose, Unit, sep = "_")])
  , 
  Induced = c(TRUE, FALSE), 
  Rep = c(1:4))

layout <- cbind(Well_locations, layout)

expanded_layout <- dcast(layout, Induced + Rep ~ Organism + Chemical, value.var = "Well")

setnames(Chemical_Data, t(Chemical_Data)[,1])

setnames(
  Chemical_Data,
  "Time [s]",
  "Well")

Chemical_Data <-
  Chemical_Data[grep("[A-z][0-9]{1,2}", Well)]

Chemical_Data <- 
  melt(
    Chemical_Data, 
    id.vars = "Well", 
    variable.name = "Time", 
    value.name = "OD600", 
    na.rm = TRUE)

Chemical_Data <- layout[Chemical_Data, on = .(Well)]

Chemical_Data[, c("Chemical", "Dose", "Unit") := tstrsplit(Chemical, "_")]

Chemical_Data[, Date := today]

Chemical_Data <- Chemical_Data[, .(Date, Well, Organism, Chemical, Dose, Unit, Induced, Rep, Time, OD600)]

Chemical_Data[, Time := as.numeric(as.character(Time))]

Chemical_Data[, Instrument := "icontrol"]

Chemical_Data <- Ryans_Chemicals[, .(Chemical, cJMP)][Chemical_Data, on = .(Chemical == Chemical)]

Chemical_Data[, OD600 := as.numeric(OD600)]

fwrite(Chemical_Data, 
       paste0("Sheets/experiment", "-", today, ".tsv"),
       sep = "\t")

Experiments <- Chemical_Data[Experiments, on = Experiments.pk]

Experiments[,
  c("Organism",
    "Chemical",
    "Dose",
    "Unit",
    "Induced",
    "Rep",
    "OD600",
    "cJMP",
    "Instrument") := .(
      i.Organism,
      i.Chemical,
      i.Dose,
      i.Unit,
      i.Induced,
      i.Rep,
      i.OD600,
      i.cJMP,
      i.Instrument)]

Experiments <- Experiments[, .SD, .SDcols = !patterns("^i.")]

dbWriteTable(chem_gen_db, "Experiments", Chemical_Data, append = TRUE)

dbDisconnect(chem_gen_db)
