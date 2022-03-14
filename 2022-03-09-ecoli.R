today = "2022-03-09"

Chemical_Data <- fread(
  paste0(today, ".tsv"), 
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
  NA,
  1, 
  30, 
  1, 
  30, 
  750)

Units <- c(
  NA,
  "%w/v",
  "ng/mL",
  "mM",
  "ug/mL",
  "mM")

source("icontrol_analyzer.R")

