require("pacman")
p_load(data.table, ggplot2, ggtext)

# Extract organism name from file path and format
get_organism <- function(file_path) {
  org <- sub(".*promoters_fluorescence_(.*)\\.tsv$", "\\1", file_path)
  formatted_org <- switch(org,
    eco = "*E. coli*",
    ecl = "*E. cloacae*",
    kpn = "*K. pneumoniae*",
    org # Default to the original if not matched
  )
  return(formatted_org)
}

# Process data from file
process_data <- function(file_path) {
  organism_name <- get_organism(file_path) # Get formatted organism name
  promoters <- fread(file_path, na.strings = "null")
  promoters[, GFP_rel := GFP / OD600]

  # Reshape data
  promoters <- dcast(promoters,
    Promoter_dCas9 + BioRep ~ Reporter + Guide + IPTG,
    value.var = "GFP_rel"
  )

  # Calculate relative GFP levels
  GFP <- promoters[, .(Promoter_dCas9,
    BioRep,
    uninduced = (GFP_gmc6_0 - empty_empty_0) / (GFP_empty_0 - empty_empty_0),
    induced = (GFP_gmc6_0.001 - empty_empty_0.001) / (GFP_empty_0.001 - empty_empty_0.001)
  )]

  # Go through the motions here, to get a standard deviation
  GFP_no_guide <- promoters[, .(Promoter_dCas9,
    BioRep,
    uninduced = GFP_empty_0 - empty_empty_0,
    induced = GFP_empty_0.001 - empty_empty_0.001
  )]


  GFP_median <- GFP[, .(
    uninduced = median(uninduced),
    induced = median(induced)
  ),
  by = .(Promoter_dCas9)
  ]

  GFP_no_guide_median <- GFP_no_guide[, .(
    uninduced = median(uninduced),
    induced = median(induced)
  )]

  # Calculate standard error of the mean (SEM)
  GFP_sem <- GFP[, .(
    uninduced = sd(uninduced) / sqrt(.N),
    induced = sd(induced) / sqrt(.N)
  ),
  by = .(Promoter_dCas9)
  ]

  GFP_no_guide_sem <- GFP_no_guide[, .(
    uninduced = sd(uninduced) / sqrt(.N),
    induced = sd(induced) / sqrt(.N)
  )]

  GFP_median <- melt(GFP_median, id.vars = "Promoter_dCas9", variable.name = "IPTG", value.name = "GFP")
  GFP_sem <- melt(GFP_sem, id.vars = "Promoter_dCas9", variable.name = "IPTG", value.name = "SEM")

  # Combine median and SEM data
  GFP_experiment <- GFP_median[GFP_sem, on = .(Promoter_dCas9, IPTG)]
  GFP_experiment[, Organism := organism_name] # Add formatted organism column

  GFP_no_guide_median[, Promoter_dCas9 := "No guide"]
  GFP_no_guide_sem[, Promoter_dCas9 := "No guide"]

  GFP_no_guide_median <- melt(GFP_no_guide_median, id.vars = "Promoter_dCas9", variable.name = "IPTG", value.name = "GFP")
  GFP_no_guide_sem <- melt(GFP_no_guide_sem, id.vars = "Promoter_dCas9", variable.name = "IPTG", value.name = "SEM")


  # Combine median and SEM data
  GFP_no_guide_experiment <- GFP_no_guide_median[GFP_no_guide_sem, on = .(Promoter_dCas9, IPTG)]
  GFP_no_guide_experiment[, Organism := organism_name] # Add formatted organism column

  GFP_experiment[, Guide := "With Guide"]
  GFP_no_guide_experiment[, Guide := "No Guide"]

  GFP_no_guide_experiment[, SEM := SEM / GFP]
  GFP_no_guide_experiment[, GFP := GFP / GFP]

  GFP_experiment <- rbind(GFP_experiment, GFP_no_guide_experiment)


  GFP_no_guide_median

  return(GFP_experiment)
}

# Process all organism files
files <- list.files("promoter_testing", pattern = "promoters_fluorescence_.*\\.tsv$", full.names = TRUE)
data <- rbindlist(lapply(files, process_data))

# Organize data for plotting
data[, Promoter_dCas9 := factor(Promoter_dCas9, levels = c("No guide", "PLlacO1", setdiff(unique(Promoter_dCas9), c("No guide", "PLlacO1"))))] # Order promoters
setorder(data, Organism, Promoter_dCas9) # Order data for plotting
data[, Promoter_Organism := paste(Promoter_dCas9, Organism, sep = "_")]


# Map X-axis positions explicitly
data[, Promoter_Organism := factor(Promoter_Organism, levels = unique(Promoter_Organism))]
data[, x_pos := as.numeric(Promoter_Organism)] # Numeric X-axis positions

ratio_data <- dcast(data, Promoter_dCas9 + Organism ~ IPTG, value.var = "GFP")
ratio_data[, ratio := uninduced / induced]


data <- merge(data, ratio_data[, .(Promoter_dCas9, Organism, ratio)], by = c("Promoter_dCas9", "Organism"), all.x = TRUE)
data[, ratio_formatted := as.character(signif(ratio, 2))]
data[ratio > 1 & IPTG == "induced", ratio_formatted := paste0("**", ratio_formatted, " x", "**")]
data[ratio <= 1 | IPTG == "uninduced", ratio_formatted := ""]

# Introduce spacing between groups
spacing <- 1 # Adjust this value to control the gap size
organism_positions <- data[, .(x_start = min(x_pos), x_end = max(x_pos)), by = Organism]
organism_positions[, x_start := x_start + (group_id <- .I - 1) * spacing] # Offset x_start for gaps
organism_positions[, x_end := x_end + (.I - 1) * spacing] # Offset x_end for gaps
data[, x_pos := x_pos + (as.integer(factor(Organism)) - 1) * spacing] # Adjust data x_pos to match


# Verify the structure of the data and grouping
print("Combined Data:")
print(head(data))
print("Organism Positions:")
print(organism_positions)

# Plot data
p <- ggplot(data, aes(x = x_pos, y = GFP, fill = IPTG, group = IPTG)) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    color = "black"
  ) +
  # Use SEM for error bars
  geom_errorbar(
    aes(ymin = GFP - SEM, ymax = GFP + SEM),
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black",
    width = 0.5
  ) +
  # Add labels for "induced" bars with the ratio values
  geom_richtext(
    data = data,
    aes(label = ratio_formatted, y = (GFP + SEM) * (1 + data[, median(GFP - SEM) * 0.25])),
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black",
    fill = NA,
    label.colour = NA,
    hjust = 0,
    # size = 6,
    angle = 90
  ) +
  # Horizontal lines under organism groups with visible breaks
  geom_segment(
    data = organism_positions,
    aes(
      x = x_start - 0.5,
      xend = x_end + 0.5,
      y = data[, min(GFP - SEM) * 0.9],
      yend = data[, min(GFP - SEM) * 0.9]
    ),
    inherit.aes = FALSE, color = "black", linewidth = 2
  ) +
  # Organism labels below the line using ggtext
  geom_richtext(
    data = organism_positions,
    fill = NA,
    label.colour = NA,
    aes(x = (x_start + x_end) / 2, y = data[, min(GFP - SEM) * 0.7], label = Organism),
    inherit.aes = FALSE,
    size = 7,
    vjust = 0.5
  ) +
  scale_x_continuous(breaks = data$x_pos, labels = data$Promoter_dCas9) +
  scale_fill_manual(values = c("grey", "#db281c")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  # write
  theme_minimal() +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = 16, color = "black"),
    plot.title = element_markdown(size = 16), # Enable markdown for the title

    # change title.x sizes but not using markdown
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),

    # legend.title = element_markdown(), # Enable markdown for legend title
    plot.caption = element_markdown() # Enable markdown for caption
  ) +
  labs(
    title = "dCas9 Promoter Efficacy on GFP Knockdown",
    x = "Promoters (Grouped by Organism)",
    y = "Relative GFP Fluorescence",
    fill = "CRISPRi Induction"
  ) +
  # capitalize the key in legends,
  guides(fill = guide_legend(title.theme = element_text(face = "bold"))) +
  # Scale continuous on log10 y
  scale_y_log10()





print(p)
