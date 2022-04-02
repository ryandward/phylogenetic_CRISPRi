---
title: "Effect of Ampicillin in 3 CRISPRi Libraries in EZ media. "
output:
  html_document: 
    toc: yes
    keep_md: yes
  html_notebook: default
  pdf_document:
    fig_height: 4
  word_document: default
editor_options: 
  chunk_output_type: console
---

# Raw Values


```r
library(pacman)
p_load(data.table, ggplot2, RSQLite)

chem_gen_db <- dbConnect(RSQLite::SQLite(), "../chem_gen.db")
Experiments <- data.table(dbReadTable(chem_gen_db, "Experiments"))
Ryan_Strains <- data.table(dbReadTable(chem_gen_db, "Ryan_Strains"))
dbDisconnect(chem_gen_db)

today <- "2022-03-19"

Experiments <- Ryan_Strains[Experiments[Date == today], on = .(Organism)]

for (i in unique(Experiments[!is.na(Species), Species])) {
  
  this.plot <- ggplot(
    Experiments[Species == i, 
                .(OD600, `Strain and Drug` = paste(Remark, ":", Chemical, Dose, Unit)), 
                by = .(Time)],
    aes(x = Time, y = OD600, color = `Strain and Drug`, fill = `Strain and Drug`)) +
    ylim(0 - 0.15, max(Experiments$OD600) + 0.15) +
    geom_smooth(formula = y ~ s(x, bs = "cs"), method = "gam") +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    ggtitle(paste(i, "growth curves on", today))
  
  plot(this.plot)}
```

![](2022-03-19-ampicillin_files/figure-html/unnamed-chunk-1-1.png)<!-- -->![](2022-03-19-ampicillin_files/figure-html/unnamed-chunk-1-2.png)<!-- -->![](2022-03-19-ampicillin_files/figure-html/unnamed-chunk-1-3.png)<!-- -->

\pagebreak

# Fitted Values using GrowthcurveR


```r
library(pacman)
p_load(data.table, ggplot2, RSQLite)

chem_gen_db <- dbConnect(RSQLite::SQLite(), "../chem_gen.db")
Fitted_Experiments <- data.table(dbReadTable(chem_gen_db, "Fitted_Experiments"))
Ryan_Strains <- data.table(dbReadTable(chem_gen_db, "Ryan_Strains"))
dbDisconnect(chem_gen_db)

today <- "2022-03-19"

Fitted_Experiments <- Ryan_Strains[Fitted_Experiments[Date == today], on = .(Organism)]

for (i in unique(Fitted_Experiments[!is.na(Species), Species])) {
  
  this.plot <- ggplot(
    Fitted_Experiments[Species == i, 
                .(OD600_fit, `Strain and Drug` = paste(Remark, ":", Chemical, Dose, Unit)), 
                by = .(Hour)],
    aes(x = Hour, y = OD600_fit, color = `Strain and Drug`, fill = `Strain and Drug`)) +
    ylim(0 - 0.15, max(Fitted_Experiments$OD600_fit) + 0.15) +
    geom_smooth(formula = y ~ s(x, bs = "cs"), method = "gam") +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    ggtitle(paste(i, "growth curves on", today))
  
  plot(this.plot)}
```

![](2022-03-19-ampicillin_files/figure-html/unnamed-chunk-2-1.png)<!-- -->![](2022-03-19-ampicillin_files/figure-html/unnamed-chunk-2-2.png)<!-- -->![](2022-03-19-ampicillin_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

\pagebreak