library(tidyverse) # Load tidyverse for data manipulation and visualization
library(data.table) # Load data.table for fast and efficient data manipulation
library(ggtext) # Load markdown for formatting text (e.g., italicizing species names)

# ============================
# Functions Section
# ============================

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
