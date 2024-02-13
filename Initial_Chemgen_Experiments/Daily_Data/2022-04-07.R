library(pacman)
p_load(data.table, RSQLite)


chem_gen_db <- dbConnect(RSQLite::SQLite(), "chem_gen.db")

Experiments <- data.table(dbReadTable(chem_gen_db, "Experiments"))
Experiments.pk <- c("Date", "Well", "Time")

Experiments[, Organism := factor(Organism)]
# Experiments[, Induced := as.logical(Induced)]

Ryans_Chemicals <- data.table(dbReadTable(chem_gen_db, "Ryans_Chemicals"))

today = "2022-04-07"

Chemical_Data <- fread(
  paste0("Daily_Data/", today, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

Organisms <- c(
  "5147",
  "5142",
  "5123",
  "5144",
  "5133",
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
    0,
    0,
    0)

Induced <-
  rep(Induced, each = 2)


Reps <- rbind(
  CJ(Well_row = rep(c(1,3,5,7,9,11),1), Rep = rep(1:8)), 
  CJ(Well_row = rep(c(2,4,6,8,10,12),1), Rep = rep(9:16)))

Well_col <- rep(c(c("A","B","C","D"),c("E","F","G","H")), 12)

Reps <- cbind(Reps, Well_col)

Reps <- Reps[, .(Well_col, Well_row, Rep)]

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


# setnames(Chemical_Data, t(Chemical_Data)[,2])

setnames(
  Chemical_Data,
  "Time [s]",
  "Well")

Chemical_Data <- Chemical_Data[-1,]

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

Chemical_Data[, Date := today]

Chemical_Data <- Chemical_Data[, .(Date, Well, Organism, Chemical, Dose, Induced, Rep, Time, OD600, Media)]

Chemical_Data[, Time := as.numeric(as.character(Time))]

Chemical_Data[, Instrument := "icontrol"]

Chemical_Data <- Ryans_Chemicals[, .(Chemical, cJMP, Unit)][Chemical_Data, on = .(Chemical == Chemical)]

Chemical_Data[, OD600 := as.numeric(OD600)]

fwrite(Chemical_Data, 
       paste0("Sheets/experiment", "-", today, ".tsv"),
       sep = "\t")

Experiments <- Experiments[Chemical_Data, on = Experiments.pk]

Experiments[,
            c("Organism",
              "Chemical",
              "Dose",
              "Unit",
              "Induced",
              "Rep",
              "OD600",
              "cJMP",
              "Instrument",
              "Media") := .(
                i.Organism,
                i.Chemical,
                i.Dose,
                i.Unit,
                i.Induced,
                i.Rep,
                i.OD600,
                i.cJMP,
                i.Instrument,
                i.Media)]

Experiments <- Experiments[, .SD, .SDcols = !patterns("^i.")]

dbWriteTable(chem_gen_db, "Experiments", Experiments, append = TRUE)

dbDisconnect(chem_gen_db)
