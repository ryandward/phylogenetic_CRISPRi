chem_gen <- fread('chem_gen.tsv')

chem_gen[, `Antibiotic/Condition` := tolower(`Antibiotic/Condition`)]

chem_gen[`Antibiotic/Condition` %like% "nickel", `Antibiotic/Condition` := "nickel(ii)"]

OMP_OMLP <- most_phenotypes_by_group_loc[STEPdb_loc_code %in% c("OMP", "OMLP")]

OMP_OMLP <- OMP_OMLP[, .(Condition = unique(Condition))]

setorder(OMP_OMLP, Condition)

OMP_OMLP[, Condition := gsub(' [\\[\\{]', ";", Condition, perl = TRUE)]
OMP_OMLP[, Condition := gsub('[\\]\\}]', ";", Condition, perl = TRUE)]
OMP_OMLP[, Condition := gsub(';;', ';', Condition)]
OMP_OMLP[, Condition := gsub(';$', '', Condition, perl = TRUE)]
OMP_OMLP[, c("Condition", "Dose", "Batch") := tstrsplit(Condition, ';', perl = TRUE)]
OMP_OMLP[, Condition := tolower(Condition)]

OMP_OMLP <- chem_gen[OMP_OMLP, on = .(`Antibiotic/Condition` == Condition)]

setnames(
  OMP_OMLP,
  c("Dose", "Antibiotic/Condition"),
  c("Shriver_Dose", "Common_Name"))

OMP_OMLP[, Batch := NULL]

fwrite(
  OMP_OMLP, 
  'shriver_omp_omlp.tsv',
  sep = "\t")

# current_chemicals <- fread('current_chemicals.tsv')
# current_chemicals[, Name := tolower(Name)]
# current_chemicals[, `alternative name` := tolower(Name)]
# 
# OMP_OMLP <- current_chemicals[, .(`cJMP###`, Name, `location (container)`, form)][OMP_OMLP[!is.na(UID)], on = .( Name == `Antibiotic/Condition`)]
# 
# setnames(
#   OMP_OMLP,
#   c("Dose", "location (container)", "Name"),
#   c("Shiver_Dose", "Location", "Common_Name"))