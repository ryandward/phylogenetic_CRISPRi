
library(tidyverse)
library(data.table)

# Read the data
data <- fread("Sequencing/Ecloacae_IMI/growth-data.tsv", header = TRUE)

# Set the column names to the first row of the data
setnames(data, as.character(unlist(data[1,])))

# Remove the first row (now that it's been used for column names)
data <- data[-1,]

# Rename the 'Time [s]' column to 'Well'
setnames(data, old = "Time [s]", new = "Well")

# Filter rows where 'Well' matches the pattern [A-z][0-9]{1,2}
data <- data[grep("[A-z][0-9]{1,2}", data$Well)]

# Reshape the data from wide to long format
data <- melt(data, id.vars = "Well", variable.name = "Seconds", value.name = "OD600", na.rm = TRUE)

# Split Well into Row and Column
data[, c("Row", "Column") := tstrsplit(Well, "(?<=^[A-Z])", perl = TRUE)]

# Convert Column to numeric
data[, Column := as.numeric(Column)]

# Define which wells are induced and which are not
data[Row %in% c("B", "C", "D"), Induced := TRUE]
data[Row %in% c("E", "F", "G"), Induced := FALSE]

message("Plate layout for CRISPRi induction")
data %>%
  select(Row, Column, Induced) %>%
  unique() %>%
  dcast(Row ~ Column, value.var = "Induced") %>%
  print()

# Function to set the doses for the dilution series going across the COLUMNS
column_dilution_series <- function(data, min_dose, max_dose, min_col, max_col, zero_col) {
  # Determine the direction of the dilution series
  direction <- ifelse(max_col > min_col, -1, 1)

  # Generate a sequence from 1 to the absolute difference between max_col and min_col
  cols <- seq(from = 1, to = abs(max_col - min_col) + 1)

  # Add or subtract the sequence from max_col depending on the direction of the dilution series
  cols <- max_col + direction * (cols - 1)

  # Set Imipenem for Columns from max_col to min_col, based on the dilution series
  lapply(cols, function(col) {
    data[Column == as.character(col), Imipenem := max_dose / 2^(abs(max_col - col))]
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
  min_dose = 0.125, min_col = 10,
  max_dose = 32, max_col = 2,
  zero_col = 11
)

message("Plate layout for drug dilution")
data %>%
  select(Row, Column, Imipenem) %>%
  unique() %>%
  dcast(Row ~ Column, value.var = "Imipenem") %>%
  print()

# Remove wells without cells
# data <- data[!is.na(Imipenem) & !is.na(Induced)]


# Remove 's' from Seconds and convert to numeric
data[, Seconds := as.numeric(gsub("s", "", Seconds))]


# Plot the growth curves
growth_curve_plot <- ggplot(
  data,
  aes(
    x = Seconds, y = OD600,
    group = interaction(Imipenem, Induced)
  )
) +
  scale_y_log10() +
  stat_summary(fun = mean, geom = "line", aes(color = factor(Imipenem), linetype = factor(Induced))) +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(Imipenem)), alpha = 0.4) +
  labs(title = "Enterobacter cloacae library growth", x = "Seconds", y = "Growth (log OD600)", color = "Imipenem", linetype = "Induced", fill = "Imipenem") +
  scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
  theme_minimal()

print(growth_curve_plot)

ecl_data <- data
