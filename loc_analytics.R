protein_dict <- fread("protein_dict.tsv")

protein_loc <- fread("protein_localization.tsv")

protein_loc <- protein_loc[protein_dict, on = .(protein)]
