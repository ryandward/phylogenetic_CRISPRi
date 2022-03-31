today = "2022-03-15"

Chemical_Data <- fread(
  paste0(today, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

Organisms <- c(
  "5147",
  "5142")

Chemicals <- c(
  "none",
  "mecillinam",
  "nickel(ii)",
  "novobiocin",
  "urea",
  "vancomycin")

Doses <- c(
  0,
  60, 
  5, 
  60, 
  500, 
  100)

source("sunrise_analyzer.R")

