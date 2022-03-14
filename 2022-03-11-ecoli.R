today = "2022-03-11"

Chemical_Data <- fread(
  paste0(today, ".tsv"), 
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
  NA,
  4, 
  120, 
  120, 
  100, 
  50)

Units <- c(
  NA,
  "ug/mL",
  "ng/mL",
  "ug/mL",
  "mM",
  "ug/mL")

source("icontrol_analyzer.R")