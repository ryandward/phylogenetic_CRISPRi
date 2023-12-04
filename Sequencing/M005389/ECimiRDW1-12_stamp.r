require("pacman")

p_load(data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, RColorBrewer, tidyverse, magrittr, ggpubr, ggallin)

# Load the data, first column is sample, second column is spacer, third column is count
data <- fread("Sequencing/M005389/ECimiRDW1-12.tsv.gz", header = FALSE, sep = "\t", col.names = c("sample", "spacer", "count"))

data <- data %>% filter(!(spacer %like% "\\*"))

design <- fread("Sequencing/M005389/design.tsv", header = TRUE, sep = "\t")

data <- data %>% inner_join(design, by = "sample")

targets <- fread("Sequencing/M005389/ECimiRDW1-12_found_spacers_target.tsv.gz", header = TRUE, sep = "\t") %>%
  mutate(type = case_when(mismatches == 0 ~ "perfect", mismatches == 1 ~ "mismatch", mismatches == "None" ~ "control"))

types <- targets %>%
  select(spacer, chr, target, mismatches) %>%
  unique() %>%
  mutate(type = case_when(mismatches == 0 ~ "perfect", mismatches == 1 ~ "mismatch", mismatches == "None" ~ "control"))

data <- data %>% inner_join(types, by = "spacer")

# calculate bottleneck metrics
bottleneck <- data %>%
  group_by(sample, type) %>%
  summarize(total = sum(count)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total = sum(total)) %>%
  ungroup() %>%
  mutate(fraction = total / sum(total)) %>%
  group_by(sample) %>%
  mutate(bottleneck = 1 / sum(fraction^2)) %>%
  ungroup() %>%
  select(sample, bottleneck) %>%
  unique()
