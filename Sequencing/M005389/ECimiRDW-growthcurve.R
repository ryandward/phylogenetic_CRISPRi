library(tidyverse)
library(data.table)
data <- fread("Sequencing/M005389/ECimiRDW-growth-data.tsv") %>% melt(id.vars = c("Seconds", "Temp"), variable.name = "Well", value.name = "OD600")

data[, c("Row", "Column") := tstrsplit(Well, "(?<=^[A-Z])", perl = TRUE)]

data[Row %in% c("B", "C", "D"), Induced := FALSE]

data[Row %in% c("E", "F", "G"), Induced := TRUE]

# Set Imipenem for Columns 11 to 3
for (col in 11:3) {
  data[Column == as.character(col), Imipenem := 32 / 2^(11 - col)]
}

# Set Imipenem for Column 2
data[Column == "2", Imipenem := 0]


data <- data[!is.na(Imipenem) & !is.na(Induced)]

# Remove 's' from Seconds and convert to numeric
data[, Seconds := as.numeric(gsub("s", "", Seconds))]

ggplot(data, aes(x = Seconds, y = log10(OD600), group = interaction(Imipenem, Induced))) +
  stat_summary(fun = mean, geom = "line", aes(color = factor(Imipenem), linetype = factor(Induced))) +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(Imipenem)), alpha = 0.2) +
  labs(title = "E. coli library growth", x = "Seconds", y = "Growth (log OD600)", color = "Imipenem", linetype = "Induced", fill = "Imipenem") +
  scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
  theme_minimal()
