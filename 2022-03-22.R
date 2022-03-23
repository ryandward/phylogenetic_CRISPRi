today = "2022-03-22"

Chemical_Data <- fread(
  paste0(today, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

Organisms <- c(
  "0",
  "0",
  "5123",
  "5123",
  "5144",
  "5144",
  "0",
  "0",
  "5123",
  "5123",
  "5144",
  "5144")



Chemicals <- c(
  "LB",
  "EZ")

Chemicals <-
  rep(Chemicals, each = 6)



Doses <- c(
  0,
  0)

Doses <-
  rep(Doses, each = 6)



Induced <-
  c(
    FALSE,
    FALSE)

Induced <-
  rep(Induced, each = 4)


Reps <- 
  rbind(
    CJ(
      Well_row = c(1,3,5,7,9,11), 
      Rep = rep(1:8)),
    CJ(
      Well_row = c(2,4,6,8,10,12),
      Rep = rep(9:16)))

Organism_key <- data.table(
  Organism = factor(Organisms, levels = unique(Organisms)))

Chemical_key <- data.table(
  Chemical = factor(Chemicals, levels = unique(Chemicals)))

Dose_key <- data.table(
  Dose = factor(Doses, levels = unique(Doses)))

Induced_key <- data.table(
  Induced = factor(Induced, levels = unique(Induced)))

Rep_key <- Reps

Chemical_key <- cbind(Chemical_key, Dose_key)

################################################################################

Well_locations <-
  CJ(
    Well_row = c(1:12), 
    Well_col = toupper(letters[1:8]))

Organism_key[, Well_row := c(1:12)]

Chemical_key[, Well_row := c(1:12)]

Induced_key[, Well_col := toupper(letters[1:8])]

Rep_key[, Well_col := rep(toupper(letters[1:8]), 12)]

layout <-
  Induced_key[
    Chemical_key[
      Organism_key[
        Rep_key[
          Well_locations, on = .(Well_row, Well_col)], 
        on = .(Well_row)], 
      on = .(Well_row)], 
    on = .(Well_col)]

layout[, Well := paste0(Well_col, Well_row)]

################################################################################

Chemical_Data <- Chemical_Data[`Raw data` %like% "[0-9]+s"]

Chemical_Data <- Chemical_Data[, -2]

Chemical_Data[, `Raw data` := gsub("s", "", `Raw data`)]

Chemical_Data <- t(Chemical_Data)

colnames(Chemical_Data) <- Chemical_Data[1, ]

Chemical_Data <- Chemical_Data[-1, ]



transposed_wells <- 
  CJ(Well_col = toupper(letters[1:8]),
     Well_row = c(1:12))[
       , .(Well = paste0(Well_col, Well_row))]


rownames(Chemical_Data) <- transposed_wells$Well

Chemical_Data <- data.table(Chemical_Data, keep.rownames = "Well")

Chemical_Data <- 
  melt(
    Chemical_Data, 
    id.vars = "Well", 
    variable.name = "Time", 
    value.name = "OD600", 
    na.rm = TRUE)

################################################################################


Chemical_Data <- layout[Chemical_Data, on = .(Well)]

Chemical_Data[, Date := today]

Chemical_Data <- Chemical_Data[, .(Date, Well, Organism, Chemical, Dose, Induced, Rep, Time, OD600)]

Chemical_Data[, Time := as.numeric(as.character(Time))]

Chemical_Data[, Instrument := "sunrise"]


################################################################################

chem_gen_db <- dbConnect(RSQLite::SQLite(), "chem_gen.db")

Experiments <- data.table(dbReadTable(chem_gen_db, "Experiments"))
Experiments.pk <- c("Date", "Well", "Time")

Experiments[, Organism := factor(Organism)]

Ryans_Chemicals <- data.table(dbReadTable(chem_gen_db, "Ryans_Chemicals"))

################################################################################

Chemical_Data <- Ryans_Chemicals[, .(Chemical, cJMP, Unit)][Chemical_Data, on = .(Chemical == Chemical)]

Chemical_Data[, OD600 := as.numeric(OD600)]

fwrite(Chemical_Data, 
       paste0("Sheets/experiment", "-", today, ".tsv"),
       sep = "\t")

dbWriteTable(chem_gen_db, "Experiments", Chemical_Data, append = TRUE)

dbDisconnect(chem_gen_db)

