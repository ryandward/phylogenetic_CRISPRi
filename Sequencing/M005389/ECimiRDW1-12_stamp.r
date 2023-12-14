require("pacman")

p_load(data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, RColorBrewer, tidyverse, magrittr, ggpubr, ggallin)

# Load the data, first column is sample, second column is spacer, third column is count
bottleneck_data <- fread(
  "Sequencing/M005389/ECimiRDW1-12.tsv.gz",
  header = FALSE, sep = "\t", col.names = c("sample", "spacer", "count")
) %>% filter(!(spacer %like% "\\*"))

bottleneck_data <- bottleneck_data %>% rbind(fread(
  "Sequencing/Freezer_Stocks/ECO_freezer.tsv.gz",
  header = TRUE, sep = "\t", col.names = c("sample", "spacer", "count")
) %>% filter(!(spacer %like% "\\*")))

bottleneck_design <- fread("Sequencing/M005389/design.tsv", header = TRUE, sep = "\t") %>%
  rbind(data.table(sample = "ECO_freezer", doublings = 0), fill = TRUE) %>%
  mutate(induction = ifelse(induced == TRUE, "induced", "uninduced")) %>%
  mutate(induction = factor(induction, levels = c("uninduced", "induced"))) %>%
  mutate(condition = paste(induction, imipenem)) %>%
  mutate(condition = ifelse(condition == "NA NA", "freezer stock", condition)) %>%
  mutate(condition = factor(condition, levels = unique(condition)))

bottleneck_data <- bottleneck_data %>% inner_join(bottleneck_design, by = "sample")

bottleneck_targets <- fread("Organisms/A_E_coli.targets.tsv.gz", header = TRUE, sep = "\t") %>%
  mutate(type = case_when(mismatches == 0 ~ "perfect", mismatches == 1 ~ "mismatch", mismatches == "None" ~ "control"))

bottleneck_essentials <- bottleneck_targets %>%
  filter(type == "mismatch") %>%
  select(locus_tag) %>%
  unique()

bottleneck_targets[locus_tag %in% bottleneck_essentials$locus_tag & type == "perfect", type := "perfect essential"]
bottleneck_targets[spacer %in% bottleneck_targets[type == "perfect essential"]$spacer & type == "perfect", type := "perfect essential"]

bottleneck_types <- bottleneck_targets %>%
  select(spacer, type) %>%
  unique()

bottleneck_data <- bottleneck_data %>%
  left_join(bottleneck_types, by = "spacer") %>%
  select(sample, spacer, count, imipenem, induced, replicate, doublings, type, condition)

bottleneck_data %>%
  dcast(spacer ~ sample, value.var = "count", fill = 0) %>%
  melt(id.vars = "spacer", value.name = "count", variable.name = "sample") %>%
  left_join(bottleneck_data %>% select(sample, imipenem, replicate, doublings, condition) %>% unique()) %>%
  left_join(bottleneck_data %>% select(spacer, type) %>% unique())

# Calculate frequencies of alleles at T0
botneck.t0 <- bottleneck_data %>%
  filter(doublings == 0) %>%
  group_by(type) %>%
  mutate(
    fi0 = count / sum(count),
    count0 = count
  ) %>%
  select(type, spacer, fi0, count0) %>%
  nest() %>%
  rename(data0 = data) %>%
  mutate(s0 = map_dbl(data0, ~ sum(.$count0)))

# Calculate frequencies of alleles at other times and compute Nb using the formula
botneck <- bottleneck_data %>%
  filter(doublings != 0) %>%
  group_by(type, condition, sample, doublings) %>%
  mutate(
    fis = count / sum(count)
  ) %>%
  nest() %>%
  mutate(
    ss = map_dbl(data, ~ sum(.$count))
  ) %>%
  full_join(botneck.t0) %>%
  mutate(
    data = map2(data, data0, inner_join)
  )

# Summary statistics
botneck <- botneck %>%
  mutate(
    data = map(data, ~ .x %>% mutate(ratio = ((fis - fi0)^2) / (fi0 * (1 - fi0)^2)))
  )

botneck <- botneck %>%
  mutate(
    f_hat = map_dbl(data, ~ sum(.$ratio)) * (1 / map_dbl(data, ~ n_distinct(.$spacer))),
    Nb = doublings / (f_hat - 1 / s0 - 1 / ss)
  )

botneck.stats <- botneck %>%
  group_by(condition, type) %>%
  summarise(
    Nb.med = median(Nb),
    Nb.range = max(Nb) - min(Nb),
    Nb.mean = mean(Nb),
    Nb.sd = sd(Nb)
  )

botneck <- botneck %>%
  mutate(
    type = factor(
      type,
      levels = c(
        "control",
        "perfect",
        "mismatch",
        "perfect essential"
      )
    )
  )

botneck.stats.plot <- botneck %>%
  ggplot(aes(x = type, y = Nb)) +
  # Draw columns for the median Nb
  stat_summary(
    fun.data = mean_se, # Automatically calculate mean and standard error
    geom = "col", # Draw error bars
    aes(fill = type)
  ) +
  stat_summary(
    fun.data = mean_se, # Automatically calculate mean and standard error
    geom = "errorbar", # Draw error bars
    width = 0.2, # Set the width of the error bars
    color = "black"
  ) +
  facet_wrap(~condition) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(1, 10^seq(0, 8)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black")) +
  labs(y = bquote(italic(N)[b]), x = "Type") +
  scale_fill_manual(
    values = c(
      "control" = "#bbbbbb",
      "perfect" = "#1F78B4",
      "mismatch" = "#FF7F00",
      "perfect essential" = "#E31A1C")) +
  guides(fill = "none")

print(botneck.stats.plot)
