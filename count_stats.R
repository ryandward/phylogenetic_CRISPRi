# Load libraries
library(pacman)
p_load(tidyverse, data.table, stringr)

# Read and preprocess data
freezer_counts <- list(
  ECO = fread("Batch1/ECO_freezer.tsv", col.names = c("spacer", "count")) %>% mutate(organism = "E. coli"),
  ECL = fread("Batch1/ECL_freezer.tsv", col.names = c("spacer", "count")) %>% mutate(organism = "E. cloacae"),
  KPN = fread("Batch1/KPN_freezer.tsv", col.names = c("spacer", "count")) %>% mutate(organism = "K. pneumoniae")) %>% 
  bind_rows() %>%
  mutate(
    discovered = grepl("\\*", spacer),
    spacer = gsub("\\*", "", spacer)
  )

# Create plot
freezer_counts %>%
  ggplot(aes(x = count)) + 
  geom_density() +
  scale_x_log10() +
  geom_label(
    data = freezer_counts %>% 
      group_by(discovered, organism) %>% 
      summarise(n = n(), total_count = sum(count), .groups = 'drop'),
    aes(
      x = 10, 
      y = Inf, 
      label = paste(
        "Unique: ", str_replace_all(format(n, big.mark = ","), " ", ""), 
        "\nSum: ", str_replace_all(format(total_count, big.mark = ","), " ", "")
      )
    ),
    vjust = 2, hjust = 0.5, colour = "black", fill = scales::alpha("white", 0.5), position = position_dodge(width = 1)
  ) +
  facet_grid(
    rows = vars(discovered), 
    cols = vars(organism), 
    scales = "free_y", 
    labeller = labeller(
      discovered = as_labeller(function(x) ifelse(x == "TRUE", "Discovered", "Designed")),
      organism = label_value
    )
  ) +
  theme_minimal()

