require('pacman'); 
p_load( data.table, scales, edgeR, statmod, gplots, ggplot2, corrplot, viridis, heatmaply, pheatmap, svglite, digest, ggcorrplot, corrr, colorspace, ggrepel)

# Read from file
promoters <- fread("promoters_fluorescence_eco.tsv", na.strings = "null")

promoters[, GFP_rel := GFP/OD600 ]

promoters <- dcast(promoters, 
									 Promoter_dCas9 + BioRep ~ Reporter + Guide + IPTG,
									 value.var = "GFP_rel")

GFP <- promoters[, .(Promoter_dCas9, 
										 BioRep, 
										 uninduced = (GFP_gmc6_0 - empty_empty_0)/(GFP_empty_0 - empty_empty_0), 
										 induced = (GFP_gmc6_0.001 - empty_empty_0.001)/(GFP_empty_0.001 - empty_empty_0.001))]

# Go through the motions here, to get a standard deviation
GFP_no_guide <- promoters[, .(Promoter_dCas9, 
															BioRep, 
															uninduced = GFP_empty_0 - empty_empty_0, 
															induced = GFP_empty_0.001 - empty_empty_0.001 )]

GFP_median <- GFP[, .(uninduced = median(uninduced), 
											induced = median(induced)), 
									by = .(Promoter_dCas9)]

GFP_no_guide_median <- GFP_no_guide[, .(
	Promoter_dCas9 = "No Guide", 
	uninduced = median(uninduced), 
	induced = median(induced))]

GFP_median <- melt(GFP_median, 
									 id.vars = "Promoter_dCas9", 
									 variable.name = "IPTG", 
									 value.name = "GFP")

GFP_no_guide_median <- melt(GFP_no_guide_median, 
														id.vars = "Promoter_dCas9", 
														variable.name = "IPTG", 
														value.name = "GFP")


GFP_sd <- GFP[, .(uninduced = sd(uninduced), 
									induced = sd(induced)), 
							by = .(Promoter_dCas9)]

GFP_no_guide_sd <- GFP_no_guide[, .(Promoter_dCas9 = "No Guide", 
																		uninduced = sd(uninduced), 
																		induced = sd(induced))]

GFP_sd <- melt(GFP_sd, 
							 id.vars = "Promoter_dCas9", 
							 variable.name = "IPTG", 
							 value.name = "SD")

GFP_no_guide_sd <- melt(GFP_no_guide_sd, 
												id.vars = "Promoter_dCas9", 
												variable.name = "IPTG", 
												value.name = "SD")

GFP_experiment <- GFP_median[GFP_sd, 
														 on = .(Promoter_dCas9, IPTG)]

GFP_control <- GFP_no_guide_median[GFP_no_guide_sd, 
																	 on = .(Promoter_dCas9, IPTG)]

# convert the control to the same scale (out of 1) as the individual promoters
GFP_control[, SD := SD / GFP]
GFP_control[, GFP := GFP / GFP]

GFP_experiment <- rbind(GFP_experiment, GFP_control)

GFP_experiment$Promoter_dCas9 = with(GFP_experiment, reorder(Promoter_dCas9, GFP + SD, min))

# ggplot magic??
p <- ggplot(GFP_experiment, aes(x = Promoter_dCas9, y = GFP, fill = IPTG, alpha = (sqrt(GFP)))) + 
	geom_bar(stat = "identity", 
					 color = "black", 
					 position = position_dodge()) +
	geom_errorbar(aes(ymin = GFP - SD, 
										ymax = GFP + SD), 
								width = 0.25, 
								position = position_dodge(0.9)) 


bio_reps <- max(promoters$BioRep)

print(
	p + labs(title = 
					 	paste("dCas9 Promoter Efficacy on GFP Knockdown in Escherichia coli\n",
					 				"Biological Replicates = ", bio_reps, ".",
					 				sep = ""), 
					 x = "dCas9 Promoter", 
					 y = "Relative GFP Fluorescence", 
					 fill = "CRISPRi Induction") +
		
		theme_classic() +
		
		scale_fill_manual( values = c( "green", "green" )) +
		
		geom_hline(yintercept = 1, linetype = "dashed", color = "black")
)

min_FC <- GFP_experiment[IPTG == "uninduced", .(uninduced = GFP - SD), by = .(Promoter_dCas9)][GFP_experiment[IPTG == "induced", .(induced = GFP + SD), by = .(Promoter_dCas9)], on = .(Promoter_dCas9)][, .(min_FC = uninduced/induced), by = .(Promoter_dCas9)]
max_FC <- GFP_experiment[IPTG == "uninduced", .(uninduced = GFP + SD), by = .(Promoter_dCas9)][GFP_experiment[IPTG == "induced", .(induced = GFP - SD), by = .(Promoter_dCas9)], on = .(Promoter_dCas9)][, .(max_FC = uninduced/induced), by = .(Promoter_dCas9)]

experiment_FC <- min_FC[max_FC, on = .(Promoter_dCas9)]

print(experiment_FC)
# figure out how to get these levels of GFP reduction onto the map programmatically.

GFP_Table <- dcast(GFP_experiment, Promoter_dCas9 ~ IPTG, value.var = "GFP")
print(GFP_Table[, .(knockdown = 1 - induced/uninduced), by = .(Promoter_dCas9)])
GFP_Knockdown <- GFP_Table[, .(knockdown = 1 - induced/uninduced), by = .(Promoter_dCas9)]
