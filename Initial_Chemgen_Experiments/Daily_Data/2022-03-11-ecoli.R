today = "2022-03-11"

Chemical_Data <- fread(
  paste0("Daily_Data/", today, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

Organisms <- c(
  "5111",
  "5142")

Chemicals <- c(
  "none",
  "ampicillin",
  "mecillinam",
  "novobiocin",
  "urea",
  "vancomycin")

Doses <- c(
  0,
  4, 
  120, 
  120, 
  100, 
  50)

source("icontrol_analyzer.R")