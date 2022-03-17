library(pacman)

p_load(
  data.table, 
  RSQLite, 
  RColorBrewer,
  growthcurver,
  ggplot2)

chem_gen_db <- dbConnect(RSQLite::SQLite(), "chem_gen.db")

Experiments <- data.table(dbReadTable(chem_gen_db, "Experiments"))
Experiments.pk <- c("Date", "Well", "Time")

Experiments[, Organism := factor(Organism)]

Well_Stats <- unique(Experiments[, .(Date, Well, cJMP, Chemical, Dose, Instrument, Rep)])
Well_Stats.pk <- c("Date", "Well")

Well_Stats <- Well_Stats[
  Experiments[, .(auc_l = SummarizeGrowth(Time/60/60, OD600)$vals$auc_l), 
              by = Well_Stats.pk], on = Well_Stats.pk]
  
Well_Stats <- Well_Stats[
  Experiments[, .(auc_e = SummarizeGrowth(Time/60/60, OD600)$vals$auc_e), 
              by = Well_Stats.pk], on = Well_Stats.pk]

Well_Stats <- Well_Stats[
  Experiments[, .(k = SummarizeGrowth(Time/60/60, OD600)$vals$k), 
              by = Well_Stats.pk], on = Well_Stats.pk]

Well_Stats <- Well_Stats[
  Experiments[, .(n0 = SummarizeGrowth(Time/60/60, OD600)$vals$n0), 
              by = Well_Stats.pk], on = Well_Stats.pk]

Well_Stats <- Well_Stats[
  Experiments[, .(r = SummarizeGrowth(Time/60/60, OD600)$vals$r), 
              by = Well_Stats.pk], on = Well_Stats.pk]



Fitted_Experiments <- 
  Well_Stats[
    , .(
      Hour = c(1:15), 
      OD600_fit = k  / ( 1 + ( ( k - n0 ) / n0 ) * exp( -r * c(1:15) ) ) ), 
    by = Well_Stats.pk]


Fitted_Experiments <- 
  unique(Experiments[, .(Chemical, Dose, cJMP, Organism, Induced, Instrument, Rep),
                     by = Well_Stats.pk])[Fitted_Experiments, on = Well_Stats.pk]

Fitted_Experiments[
  is.na(OD600_fit), 
  OD600_fit := 0]

Fitted_Experiments[, Dose := factor(Dose)]
Fitted_Experiments[, Induced := factor(Induced)]

dbWriteTable(chem_gen_db, 
             "Well_Stats", 
             Well_Stats, 
             overwrite = TRUE)

dbWriteTable(chem_gen_db, 
             "Fitted_Experiments", 
             Fitted_Experiments, 
             overwrite = TRUE)

dbDisconnect(chem_gen_db)


for (i in Fitted_Experiments[, unique(Date)]) {
  
  this.Fitted_Plot <- 
    ggplot(
      Fitted_Experiments[
        Date == i, 
        .(OD600_fit, Chemical = paste(Chemical, Dose)), by = .(Hour, Organism, Induced)],
      aes(
        x = Hour, 
        y = OD600_fit, 
        color = Chemical, 
        fill = Chemical)) +
    geom_smooth(
      method = "gam") +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    ggtitle(paste("Growth Curves on", i))
  
  plot(this.Fitted_Plot)
  
}
