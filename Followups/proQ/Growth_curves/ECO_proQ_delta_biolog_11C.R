library(tidyverse) # Load tidyverse for data manipulation and visualization
library(data.table) # Load data.table for fast and efficient data manipulation
library(ggtext) # Load markdown for formatting text (e.g., italicizing species names)

# Read the data file and perform some initial data transformations
data <- fread("Followups/proQ/Growth_curves/ECO_proQ_delta_biolog_11C_20241125.tsv") %>%
  select(-`Cycle Nr.`, -`Temp. [°C]`) %>%
  rename(Seconds = `Time [s]`) %>%
  melt(
    id.vars = "Seconds",
    variable.name = "Well",
    value.name = "OD600"
  )

biolog_conditions <- fread("Followups/proQ/Growth_curves/Biolog_11C.tsv")

# Create a growth curve plot to visualize OD600 over time
growth_curve_plot <- data %>%
  inner_join(biolog_conditions, by = c("Well" = "Well")) %>%
  mutate(
    `Dose Level` = factor(`Dose`),
    Hours = Seconds / 60 / 60
  ) %>%
  ggplot(
    aes(
      x = Hours,
      y = OD600,
      color = `Dose Level`
    )
  ) +
  scale_y_log10() + # Log-transform the y-axis to better visualize growth differences
  geom_line() +
  facet_wrap(~Chemical, nrow = 4) +
  theme_minimal() +
  theme(
    plot.title = element_markdown() # Enable markdown in the plot title
  ) +
  ggtitle("Growth Curves of *E. coli* WT in Biolog 11C Plates")

print(growth_curve_plot)
