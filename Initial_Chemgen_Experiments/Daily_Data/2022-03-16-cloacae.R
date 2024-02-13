today = "2022-03-16"

Chemical_Data <- fread(
  paste0("Daily_Data/", today, ".tsv"), 
  header = TRUE,
  na.strings = "NA")

Organisms <- c(
  "5123",
  "5144")

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

