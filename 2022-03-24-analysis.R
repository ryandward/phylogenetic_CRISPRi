today <- "2022-03-24"

Ryan_Strains[, Organism := factor(Organism)]

Experiments <-
  Ryan_Strains[Experiments, on = .(Organism)]

this.Fitted_Plot <- 
  ggplot(
    Experiments[
      Date == today & Species == 'Escherichia coli K-12', 
      .(OD600,
        Chemical = paste(Induced)), 
      by = .(Time)],
    aes(
      x = Time, 
      y = OD600, 
      color = Chemical, 
      fill = Chemical)) +
  geom_smooth(
    method = "gam") +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  ggtitle(paste("Growth Curves on", i))

plot(this.Fitted_Plot)