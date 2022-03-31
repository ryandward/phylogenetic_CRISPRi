today = "2022-03-09"

Chemical_Data <- fread(
  paste0("Daily_Data/", today, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

Organisms <- c(
  "5111",
  "5142")

Chemicals <- c(
  "none",
  "deoxycholate",
  "mecillinam",
  "nickel(ii)",
  "novobiocin",
  "urea")

Doses <- c(
  0,
  1, 
  30, 
  1, 
  30, 
  750)

source("icontrol_analyzer.R")

