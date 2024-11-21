source("Sequencing/IMI_Ecloacae/growthcurve-2.R")

eco_data <- data %>% mutate(Organism = "E. coli")

source("Sequencing/IMI_Kpneumoniae/growthcurve2A.R")

kpn_data <- data %>% mutate(Organism = "K. pneumoniae")

source("Sequencing/IMI_Ecoli/ECimiRDW-growthcurve.R")

ecl_data <- data %>% mutate(Organism = "E. cloacae")


all_data <- rbind(eco_data, kpn_data, ecl_data, fill = TRUE) %>% filter(!is.na(Imipenem))

median_dose <- all_data %>%
  pull(Imipenem) %>%
  median()

range_summary <- function(x) {
  data.frame(
    ymin = min(x),
    ymax = max(x)
  )
}

custom_palette <- viridis(6, option = "D")[-8]


ggplot(
  all_data %>%
    filter(followup %like% "Moved") %>%
    mutate(
      Hours = Seconds / 3600,
      Organism = paste0("*", Organism, "*"),
    ),
  aes(
    x = Hours, y = OD600,
    group = interaction(Imipenem, Induced)
  )
) +
  scale_y_log10() +
  stat_summary(fun = mean, geom = "line", aes(color = factor(Imipenem), linetype = factor(Induced)), linewidth = 0.75) +
  stat_summary(fun.data = range_summary, geom = "ribbon", aes(fill = factor(Imipenem)), alpha = 0.6) +
  labs(
    x = "Time (hours)",
    y = "Growth (OD<sub>600</sub>)",
    color = "Imipenem",
    linetype = "Induced", fill = "Imipenem"
  ) +
  scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
  scale_color_discrete(guide = guide_legend(row = ncol)) +
  facet_grid(Induced ~ Organism) +
  # use viridis color palette
  scale_color_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  theme_minimal() +
  theme(
    legend.title = element_markdown(),
    plot.title = element_markdown(),
    strip.text.y = element_text(angle = 0),
    strip.text = element_markdown(size = 20),
    # change the size of the y-axis text
    axis.text.y = element_text(size = 16),
    # change the size of the x-axis text
    axis.text.x = element_text(size = 16),
    # element markdown in axis labels
    axis.title = element_markdown(size = 20),
    axis.title.y = element_markdown(size = 20)
  )
