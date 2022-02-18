library(pacman)

p_load(data.table)

#https://www.cell.com/fulltext/S0092-8674(10)01374-7
typas_data <- fread("mmc2.tsv")

typas_data <- melt(typas_data, 
                   variable.name = "condition", 
                   value = "score", 
                   id.vars = "Gene")

typas_data[, c("ECK_name", paste0(rep("trash",5), c(1:5))) := tstrsplit(Gene, "-", type.convert = TRUE, fixed = TRUE)]

typas_conversion <- melt(typas_data[, .(ECK_name, trash1, trash2, trash3, trash4, trash5)], id.vars = "ECK_name", variable.name = "trash", value.name = "synonym")[!is.na(synonym)]

typas_conversion[, trash := NULL]

typas_data[, paste0(rep("trash",5), c(1:5)) := NULL]

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