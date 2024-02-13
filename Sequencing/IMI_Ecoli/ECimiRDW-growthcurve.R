library(tidyverse)
library(data.table)

data <- fread("Sequencing/IMI_Ecoli/ECimiRDW-growth-data.tsv") %>% melt(id.vars = c("Seconds", "Temp"), variable.name = "Well", value.name = "OD600")

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

data[Row %in% c("A", "H"), Imipenem := NA_integer_]
data[Row %in% c("A", "H"), Imipenem := NA]

data <- data[!Column %in% c(1,12)]


data <- data[
  Imipenem %in% c(0, 0.125, 0.25) & !(is.na(Induced)),
  followup := "Moved to Second Competition"
]

data <- data[is.na(followup), followup := "96-well Plate Only"] 

# Function to find the largest divisor closest to the square root
find_ncol <- function(n) {
  sqrt_n <- floor(sqrt(n))
  for (i in sqrt_n:2) {
    if (n %% i == 0) {
      return(i)
    }
  }
  return(sqrt_n)
}

# Number of levels in Imipenem
n_levels <- length(data$Imipenem %>% unique())

# Find the best number of columns
ncol <- find_ncol(n_levels)

# Create the plot
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
  labs(title = "E. coli library growth", 
    x = "Seconds",
    y = "Growth (log OD600)",
    color = "Imipenem",
    linetype = "Induced", fill = "Imipenem"
  ) +
  scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
  # scale_color_discrete(guide = guide_legend(ncol = ncol)) +
  facet_wrap(followup ~ .) +
  theme_minimal()

print(growth_curve_plot)

endpoint <- data %>%
  filter(Seconds == max(Seconds)) %>%
  dcast(Row ~ Column, value.var = "OD600")

endpoint_mat <- as.matrix(endpoint[, -1])

rownames(endpoint_mat) <- endpoint$Row

# draw matrix in ggplot with dose and induction in the boxes
library(ggplot2)


print(growth_curve_plot)
