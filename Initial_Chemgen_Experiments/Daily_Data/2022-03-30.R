today = "2022-03-30"

Chemical_Data <- fread(
  paste0("Daily_Data/", today, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

Organisms <- c(
  "163",
  "5142",
  "269",
  "5144",
  "241",
  "5146")

Organisms <-
  rep(Organisms, each = 2)

Medias <- c(
  "EZ")

Medias <-
  rep(Medias, each = 12)



Doses <- c(
  0)

Doses <-
  rep(Doses, each = 12)



Induced <-
  c(
    0,
    0.01,
    0.1,
    1)

Induced <-
  rep(Induced, each = 2)


Reps <- rbind(
  CJ(Well_col = rep(c("A","C","E","G"),6), Rep = rep(1:2)), 
  CJ(Well_col = rep(c("B","D","F","H"),6), Rep = rep(3:4)))

Well_row <- rep(c(c(1,3,5,7,9,11),c(2,4,6,8,10,12)), 8)

Reps <- cbind(Reps, Well_row)

setorder(Reps, Well_col, Well_row)

Organism_key <- data.table(
  Organism = factor(Organisms, levels = unique(Organisms)))

Media_key <- data.table(
  Media = factor(Medias, levels = unique(Medias)))

Dose_key <- data.table(
  Dose = factor(Doses, levels = unique(Doses)))

Induced_key <- data.table(
  Induced = factor(Induced, levels = unique(Induced)))

Rep_key <- Reps

Media_key <- cbind(Media_key, Dose_key)

################################################################################

Well_locations <-
  CJ(
    Well_row = c(1:12), 
    Well_col = toupper(letters[1:8]))

Organism_key[, Well_row := c(1:12)]

Media_key[, Well_row := c(1:12)]

Induced_key[, Well_col := toupper(letters[1:8])]

layout <-
  Induced_key[
    Media_key[
      Organism_key[
        Rep_key, 
        on = .(Well_row)], 
      on = .(Well_row)], 
    on = .(Well_col)]

layout[, Chemical := "none"]

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

Chemical_Data <- Chemical_Data[, .(Date, Well, Organism, Media, Dose, Induced, Rep, Time, OD600)]

Chemical_Data[, Time := as.numeric(as.character(Time))]

Chemical_Data[, Instrument := "sunrise"]

Chemical_Data[, Chemical := "none"]

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

