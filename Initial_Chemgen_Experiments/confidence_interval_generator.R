source("summarize_well_stats.R")

p_load(pracma)

target_confidence <- 0.95

Replicate_Stats <- 
  Fitted_Experiments[
    , .(OD600_fit), 
    by = .(Chemical, Dose, Organism, Date, Hour, Induced, Rep, Well)]

Grouped_Stats <-
  Replicate_Stats[
    , .(mean_OD600_fit = mean(OD600_fit),
        sd_OD600_fit = pracma::std(OD600_fit),
        sem_OD600_fit = pracma::std(OD600_fit)/sqrt(.N - 1),
        N = .N),
    by = .(Chemical, Dose, Organism, Date, Hour, Induced)]

Grouped_Stats[
  , CI_lower := mean_OD600_fit + qt((1 - target_confidence)/2, df = N - 1) * sem_OD600_fit]

Grouped_Stats[
  , CI_upper := mean_OD600_fit + qt((1 + target_confidence)/2, df = N - 1) * sem_OD600_fit]



i <- "mecillinam"
j <- "5147"
plot_shade <- "red"

plot_object <- 
  ggplot(
    Grouped_Stats[Chemical %in% c("none", i) & Induced == 1 &  Organism == j &  Date == '2022-03-15'], 
    aes(x = Hour, y = mean_OD600_fit, group = Dose, color = Dose)) +
  geom_line(
    data = Replicate_Stats[Chemical %in% c("none", i) & Induced == 1 & Organism == j &  Date == '2022-03-15'], 
    aes(x = Hour, y = OD600_fit, group = Dose), 
    color = "dark grey") +
  geom_line(
    size = 1.5, 
    alpha = 0.8, 
    color = "black") +
  geom_ribbon(
    aes(ymin = CI_lower, ymax = CI_upper, fill = Dose),
    # fill = Dose, 
    alpha = 0.2) + 
  xlim(
    Grouped_Stats[, min(Hour)], 
    Grouped_Stats[, max(Hour)]) +
  ylim(
    Grouped_Stats[, min(CI_lower, na.rm = T)], 
    Grouped_Stats[, max(CI_upper, na.rm = T)]) +
  ggtitle(i) +
  xlab("Hours") +
  ylab("Density")

print(plot_object)