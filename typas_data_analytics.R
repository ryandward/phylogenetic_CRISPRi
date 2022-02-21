library(pacman)

p_load(data.table, pheatmap)

#https://www.cell.com/fulltext/S0092-8674(10)01374-7
typas_data <- fread("mmc2.tsv")

typas_data <- melt(typas_data, 
                   variable.name = "condition", 
                   value = "score", 
                   id.vars = "Gene")

typas_data[, condition := gsub('UNSPECIFIED', "", condition)]

typas_data[, condition := gsub(' - $', "", condition)]
typas_data[, condition := gsub(' -$', "", condition)]



typas_data[, c("ECK_name", paste0(rep("extra",5), c(1:5))) := tstrsplit(Gene, "-", type.convert = TRUE, fixed = TRUE)]

typas_conversion <- melt(typas_data[, .(ECK_name, extra1, extra2, extra3, extra4, extra5)], id.vars = "ECK_name", variable.name = "extra", value.name = "synonym")[!is.na(synonym)]

typas_conversion[, extra := NULL]

typas_data[, paste0(rep("extra",5), c(1:5)) := NULL]

typas_data[, Gene := NULL]

conversion <- fread(
  "ECK to b.txt", 
  na.strings = "")[, .(long_name = V1, b_name = V3, ECK_name = V4, synonyms = V6)]

conversion[, paste0(rep("syn",19), c(1:19)) := tstrsplit(synonyms, ",")]

conversion[, synonyms := NULL]

conversion <- melt(conversion, id.vars = c("long_name", "b_name", "ECK_name"), variable.name = "synonym_number", value.name = "synonym")

conversion <- conversion[!is.na(synonym)] 

conversion[, synonym_number := NULL]

conversion <- conversion[synonym %like% "ECK"]

conversion <- melt(conversion, id.vars = c("long_name", "b_name"), variable.name = "name_type", value.name = "ECK_name")

conversion[, name_type := NULL]

conversion <- unique(conversion)

conversion_clues <- unique(typas_conversion[typas_data[!ECK_name %in% conversion$ECK_name, .(ECK_name = unique(ECK_name))], on = .(ECK_name)])

conversion_clues[!is.na(synonym), synonym := paste0(tolower(substr(synonym, 1, 3)), substr(synonym, 4, length(synonym)))]

conversion_clues[
  synonym %in% conversion$long_name, 
  b_name := conversion[long_name %in% conversion_clues$synonym][
    conversion_clues[synonym %in% conversion$long_name], 
    on = .(long_name == synonym)]$b_name]

conversions_found <- conversion_clues[!is.na(b_name)]

setnames(conversions_found, old = "synonym", new = "long_name")

conversion <- rbind(conversion, conversions_found)

typas_data <- conversion[typas_data, on = .(ECK_name)][!is.na(b_name)]

# get rid of non-unique b names
typas_data <- 
  typas_data[
    !b_name %in% typas_data[, .N, by = .(b_name)][N != nrow(typas_data[, .(unique(condition))])]$b_name]

# get rid of non-unique ECK names
typas_data <- 
  typas_data[
    !b_name %in% typas_data[, .N, by = .(ECK_name)][N != nrow(typas_data[, .(unique(condition))])]$b_name]

typas_data <- dcast(typas_data, b_name ~ condition, value.var = "score")

typas_data_matrix <- data.matrix(typas_data[, -1])

rownames(typas_data_matrix) <- typas_data$b_name

typas_data_matrix <- typas_data_matrix[complete.cases(typas_data_matrix),]

typas_cor <- cor(typas_data_matrix, use = "complete.obs")

plot_matrix <- typas_cor

break_halves <- length(unique(as.vector(plot_matrix)))

breaks <- c(
  seq(min(plot_matrix),
      median(plot_matrix),
      length.out = break_halves)[-break_halves],
  seq(median(plot_matrix),
      max(plot_matrix),
      length.out = break_halves))

breaks <- breaks[-length(breaks)]

breaks <- c(breaks, 0.99999999)

plot_colors <-
  c(colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves],
    colorRampPalette(c("white", "#007ac1"))(break_halves)[-c(1, break_halves)])

plot_colors <-
  colorRampPalette(c("white", "#007ac1"))(break_halves * 2 - 1)[-c(1, break_halves)]

plot_colors <- c(plot_colors, "dark grey")

to_plot_title <- paste("Raw Count Condition Correlations")

to_plot <- pheatmap(
  plot_matrix,
  col = plot_colors,
  breaks = breaks,
  border_color = NA,
  # cellwidth = 20,
  # cellheight = 20,
  main = to_plot_title,
  angle_col = 315,
  # fontsize_col = 10,
  # fontsize_row = 10,
  # cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  clustering_method = "ward.D2",
  # clustering_distance_rows = "maximum",
  # clustering_distance_cols = "maximum"
)

plot(to_plot$tree_row, cex = 0.35)
abline(h=10, col="red", lty=2, lwd=2)

condition_groups <- data.matrix(sort(cutree(to_plot$tree_row, h=10)))
condition_groups <- data.table(condition_groups, keep.rownames = "condition")
setnames(condition_groups, "V1", "group")

fwrite(condition_groups, "condition_groups.tsv", sep = "\t")