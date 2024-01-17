library(tidyverse)
library(data.table)

data <- fread("Sequencing/M005389/ECimiRDW-growth-data.tsv") %>% melt(id.vars = c("Seconds", "Temp"), variable.name = "Well", value.name = "OD600")

# Split Well into Row and Column
data[, c("Row", "Column") := tstrsplit(Well, "(?<=^[A-Z])", perl = TRUE)]

# Convert Column to numeric
data[, Column := as.numeric(Column)]

# Define which wells are induced and which are not
data[Row %in% c("B", "C", "D"), Induced := FALSE]
data[Row %in% c("E", "F", "G"), Induced := TRUE]

message("Plate layout for CRISPRi induction")
data %>%
  select(Row, Column, Induced) %>%
  unique() %>%
  dcast(Row ~ Column, value.var = "Induced") %>%
  print()

# Function to set the doses for the dilution series going across the COLUMNS
column_dilution_series <- function(data, min_dose, max_dose, min_col, max_col, zero_col) {
  # Set Imipenem for Columns from max_col to min_col, based on the dilution series
  lapply(max_col:min_col, function(col) {
    data[Column == as.character(col), Imipenem := max_dose / 2^(max_col - col)]
  })

  # Set the specified column to zero
  data[Column == as.character(zero_col), Imipenem := 0]

  # Check if the doses are aligned with the columns
  if (
    any(data[Column == as.character(min_col), Imipenem] != min_dose) ||
      any(data[Column == as.character(max_col), Imipenem] != max_dose) ||
      any(data[Column == as.character(zero_col), Imipenem] != 0)
  ) {
    stop("Error: The doses are not aligned with the columns.")
  }

  return(data)
}

# Call the function with explicit parameters
data <- column_dilution_series(
  data = data,
  min_dose = 0.125, min_col = 3,
  max_dose = 32, max_col = 11,
  zero_col = 2
)

message("Plate layout for drug dilution")
data %>%
  select(Row, Column, Imipenem) %>%
  unique() %>%
  dcast(Row ~ Column, value.var = "Imipenem") %>%
  print()

# Remove wells without cells
data <- data[!is.na(Imipenem) & !is.na(Induced)]


# Remove 's' from Seconds and convert to numeric
data[, Seconds := as.numeric(gsub("s", "", Seconds))]


# Plot the growth curves
growth_curve_plot <- ggplot(
  data,
  aes(
    x = Seconds, y = log10(OD600),
    group = interaction(Imipenem, Induced)
  )
) +
  stat_summary(fun = mean, geom = "line", aes(color = factor(Imipenem), linetype = factor(Induced))) +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(Imipenem)), alpha = 0.2) +
  labs(title = "E. coli library growth", x = "Seconds", y = "Growth (log OD600)", color = "Imipenem", linetype = "Induced", fill = "Imipenem") +
  scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
  theme_minimal()

print(growth_curve_plot)
