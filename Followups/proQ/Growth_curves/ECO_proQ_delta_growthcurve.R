library(tidyverse) # Load tidyverse for data manipulation and visualization
library(data.table) # Load data.table for fast and efficient data manipulation
library(ggtext) # Load markdown for formatting text (e.g., italicizing species names)

# ============================
# Functions Section
# ============================

# Function to set the doses for the dilution series across the columns of the plate
column_dilution_series <- function(data, min_dose, max_dose, min_col, max_col, zero_col) {
  # Determine the direction of the dilution series (increasing or decreasing)
  direction <- ifelse(max_col > min_col, -1, 1)

  # Generate a sequence of column numbers between min_col and max_col
  cols <- seq(from = 1, to = abs(max_col - min_col) + 1)

  # Adjust the column sequence based on the direction of dilution
  cols <- max_col + direction * (cols - 1)

  # Assign Imipenem doses for each column in the dilution series
  lapply(cols, function(col) {
    data[Column == as.character(col), Imipenem := max_dose / 2^(abs(max_col - col))]
  })

  # Set the specified column's dose to zero (control column)
  data[Column == as.character(zero_col), Imipenem := 0]

  # Check if the doses are aligned correctly
  if (
    any(data[Column == as.character(min_col), Imipenem] != min_dose) ||
      any(data[Column == as.character(max_col), Imipenem] != max_dose) ||
      any(data[Column == as.character(zero_col), Imipenem] != 0)
  ) {
    stop("Error: The doses are not aligned with the columns.")
  }

  return(data)
}

# Function to find the largest divisor closest to the square root of n
find_ncol <- function(n) {
  sqrt_n <- floor(sqrt(n))
  for (i in sqrt_n:2) {
    if (n %% i == 0) {
      return(i)
    }
  }
  return(sqrt_n)
}

# ============================
# Main Code Section
# ============================

# Read the data file and perform some initial data transformations
data <- fread("Followups/proQ/Growth_curves/proQ_delta_growth.tsv") %>%
  select(-`Cycle Nr.`, -`Time [s]`) %>% # Remove unwanted columns for simplicity
  melt(id.vars = c("Time (hours)", "Temp. [°C]"), variable.name = "Well", value.name = "OD600") %>%
  # Rename columns for easier handling
  rename(Hours = `Time (hours)`, Temp = `Temp. [°C]`)

# Split the "Well" column into separate Row and Column identifiers
data[, c("Row", "Column") := tstrsplit(Well, "(?<=^[A-Z])", perl = TRUE)]

# Convert the "Column" from character to numeric for easy numerical operations
data[, Column := as.numeric(Column)]

# Label specific rows as either "WT" (wild type) or "ΔproQ" (mutant) based on their position
data[Row %in% c("B", "C", "D"), Strain := "WT"] # Rows B, C, D are WT
data[Row %in% c("E", "F", "G"), Strain := "Δ*proQ*"] # Rows E, F, G are ΔproQ

# Convert the "Strain" column to a factor to explicitly define it as categorical
data[, Strain := as.factor(Strain)]

# Display a message with a preview of the plate layout for CRISPRi induction
message("Plate layout")
data %>%
  select(Row, Column, Strain) %>%
  unique() %>%
  dcast(Row ~ Column, value.var = "Strain") %>% # Reshape data to show the plate layout
  print()

# Apply the dilution series function to the data with specified parameters
data <- column_dilution_series(
  data = data,
  min_dose = 0.125, # Minimum Imipenem concentration
  min_col = 10, # Column where minimum concentration is applied
  max_dose = 32, # Maximum Imipenem concentration
  max_col = 2, # Column where maximum concentration is applied
  zero_col = 11 # Control column with zero Imipenem
)

# Display the plate layout with the Imipenem doses
message("Plate layout for drug dilution")
data %>%
  select(Row, Column, Imipenem) %>%
  unique() %>%
  dcast(Row ~ Column, value.var = "Imipenem") %>% # Reshape to visualize the plate layout
  print()

# Remove wells that do not have cells (e.g., empty or control wells)
data <- data[!is.na(Imipenem) & !is.na(Strain)]

# Remove 's' from Hours values and convert to numeric for plotting purposes
data[, Hours := as.numeric(gsub("s", "", Hours))]

# Set Imipenem to NA for rows A and H, which are assumed to be empty or irrelevant controls
data[Row %in% c("A", "H"), Imipenem := NA_integer_]
data[Row %in% c("A", "H"), Imipenem := NA]

# Remove columns 1 and 12 which are assumed to be empty or not needed for analysis
data <- data[!Column %in% c(1, 12)]

# Identify follow-up wells that should move to the "Second Competition" phase
data <- data[
  Imipenem %in% c(0, 0.125, 0.25) & !(is.na(Strain)),
  followup := "Moved to Second Competition"
]

# Label the remaining wells as "96-well Plate Only"
data <- data[is.na(followup), followup := "96-well Plate Only"]

# Calculate the number of unique Imipenem concentrations used
n_levels <- length(data$Imipenem %>% unique())

# Find the best number of columns for plotting based on the number of levels
ncol <- find_ncol(n_levels)

# Create a growth curve plot to visualize OD600 over time
growth_curve_plot <- ggplot(
  data,
  aes(
    x = Hours, y = OD600,
    group = interaction(Imipenem, Strain) # Group by unique combinations of Imipenem and Strain
  )
) +
  scale_y_log10() + # Log-transform the y-axis to better visualize growth differences
  # Plot the mean growth curve for each group
  stat_summary(fun = mean, geom = "line", aes(color = factor(Imipenem), linetype = factor(Strain))) +
  # Add ribbons to indicate standard error around the mean
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(Imipenem)), alpha = 0.4) +
  # Add labels to the plot
  labs(
    title = "*E. coli* growth in WT and Δ*proQ* strains with varying Imipenem concentrations",
    x = "Hours",
    y = "Growth (log OD600)",
    color = "Imipenem",
    linetype = "Strain",
    fill = "Imipenem"
  ) +
  # Customize line types for the different strains
  scale_linetype_manual(values = c("Δ*proQ*" = "dashed", "WT" = "solid")) +
  # Facet by strain to separate the plots for better visualization
  facet_grid(~Strain) +
  theme_minimal() + # Use a minimal theme for a clean look
  theme(
    plot.title = element_markdown(), # Enable markdown in the plot title
    axis.title.x = element_markdown(), # Enable markdown in the x-axis title
    axis.title.y = element_markdown(), # Enable markdown in the y-axis title
    strip.text = element_markdown(), # Enable markdown in the facet labels
    legend.text = element_markdown() # Enable markdown in the legend labels
  )

# Print the growth curve plot
print(growth_curve_plot)
