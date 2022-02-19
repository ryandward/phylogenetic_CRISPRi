library(pacman)

p_load(data.table, pheatmap)

#https://www.cell.com/fulltext/S0092-8674(10)01374-7
typas_data <- fread("mmc2.tsv")

typas_data <- melt(typas_data, 
                   variable.name = "condition", 
                   value = "score", 
                   id.vars = "Gene")

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

