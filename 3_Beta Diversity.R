# Pairwise comparisons-PERMANOVA-beta diversity between types, varieties and biomass
library(pairwiseAdonis)
library(vegan)
library(openxlsx)
library(phyloseq)
library(ggplot2)
library(ggpubr)

setwd("~path")
fungi_rare <- readRDS('fungi_rarefied.rds')
bacteria_rare <- readRDS('bacteria_rarefied.rds')

# Start with pairwise comparisons between varieties in fungi and yeasts for each type
# Subset samples to the desired Type
fun_C <- subset_samples(fungi_rare, Type == "C")
fun_GU <- subset_samples(fungi_rare, Type == "GU")
fun_GR <- subset_samples(fungi_rare, Type == "GR")

# Transform raw counts to relative abundances for each type
fun_C_100 <- transform_sample_counts(fun_C, function(OTU) 100*OTU/sum(OTU))
fun_GU_100 <- transform_sample_counts(fun_GU, function(OTU) 100*OTU/sum(OTU))
fun_GR_100 <- transform_sample_counts(fun_GR, function(OTU) 100*OTU/sum(OTU))

# Interactions
mycmpfactor_C_100 <- interaction(data.frame(fun_C_100@sam_data)$Variety)
mycmpfactor_GU_100 <- interaction(data.frame(fun_GU_100@sam_data)$Variety)
mycmpfactor_GR_100 <- interaction(data.frame(fun_GR_100@sam_data)$Variety)

# Run pairwise comparisons for varieties in each type
mympairwiseperm_C_100 <- pairwise.adonis(fun_C_100@otu_table, mycmpfactor_C_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
mympairwiseperm_GU_100 <- pairwise.adonis(fun_GU_100@otu_table, mycmpfactor_GU_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
mympairwiseperm_GR_100 <- pairwise.adonis(fun_GR_100@otu_table, mycmpfactor_GR_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

# Save files
write.xlsx(data.frame(mympairwiseperm_C_100), file="C:/Users/User/Desktop/Pairwise_Fungi_C.xlsx")
write.xlsx(data.frame(mympairwiseperm_GU_100), file="C:/Users/User/Desktop/Pairwise_Fungi_GU.xlsx")
write.xlsx(data.frame(mympairwiseperm_GR_100), file="C:/Users/User/Desktop/Pairwise_Fungi_GR.xlsx")

##----------------------------------------------------------------------------------------------------------
# Continue with bacteria
# Subset samples to the desired Type
bac_C <- subset_samples(bacteria_rare, Type == "C")
bac_GU <- subset_samples(bacteria_rare, Type == "GU")
bac_GR <- subset_samples(bacteria_rare, Type == "GR")

# Transform raw counts to relative abundances in each type
bac_C_100 <- transform_sample_counts(bac_C, function(OTU) 100*OTU/sum(OTU))
bac_GU_100 <- transform_sample_counts(bac_GU, function(OTU) 100*OTU/sum(OTU))
bac_GR_100 <- transform_sample_counts(bac_GR, function(OTU) 100*OTU/sum(OTU))

# Interactions
mycmpfactor_C_100b <- interaction(data.frame(bac_C_100@sam_data)$Variety)
mycmpfactor_GU_100b <- interaction(data.frame(bac_GU_100@sam_data)$Variety)
mycmpfactor_GR_100b <- interaction(data.frame(bac_GR_100@sam_data)$Variety)

# Run pairwise comparisons between varieties in each type
mympairwiseperm_C_100b <- pairwise.adonis(bac_C_100@otu_table, mycmpfactor_C_100b, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
mympairwiseperm_GU_100b <- pairwise.adonis(bac_GU_100@otu_table, mycmpfactor_GU_100b, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
mympairwiseperm_GR_100b <- pairwise.adonis(bac_GR_100@otu_table, mycmpfactor_GR_100b, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

# Save files
write.xlsx(data.frame(mympairwiseperm_C_100b), file="C:/Users/User/Desktop/Pairwise_Bacteria_C.xlsx")
write.xlsx(data.frame(mympairwiseperm_GU_100b), file="C:/Users/User/Desktop/Pairwise_Bacteria_GU.xlsx")
write.xlsx(data.frame(mympairwiseperm_GR_100b), file="C:/Users/User/Desktop/Pairwise_Bacteria_GR.xlsx")

##-------------------------------------------------------------------------------------------------------
# All types collectively, not taking account each variety separately
# Transform raw counts to relative abundance both for fungi and bacteria
all_Types_100f <- transform_sample_counts(fungi_rare, function(OTU) 100*OTU/sum(OTU))
all_Types_100b <- transform_sample_counts(bacteria_rare, function(OTU) 100*OTU/sum(OTU))

# Interaction
mycmpfactor_all_100f <- interaction(data.frame(all_Types_100f@sam_data)$Type, data.frame(all_Types_100f@sam_data)$Type)
mycmpfactor_all_100b <- interaction(data.frame(all_Types_100b@sam_data)$Type, data.frame(all_Types_100b@sam_data)$Type) 

# Pairwise comparison between Types both in fungi and bacteria
mympairwiseperm_all_100f <- pairwise.adonis(all_Types_100f@otu_table, mycmpfactor_all_100f, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
mympairwiseperm_all_100b <- pairwise.adonis(all_Types_100b@otu_table, mycmpfactor_all_100b, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

# Save xlsx files
write.xlsx(data.frame(mympairwiseperm_all_100f), file="C:/Users/User/Desktop/Pairwise_All_Types_Fungi.xlsx")
write.xlsx(data.frame(mympairwiseperm_all_100b), file="C:/Users/User/Desktop/Pairwise_All_Types_Bacteria.xlsx")

##-------------------------------------------------------------------------------------------------------
# PERMANOVA for factor = Types
mypermanova_All_100f <- adonis2(all_Types_100f@otu_table ~ Type, method = "bray", data = data.frame(all_Types_100f@sam_data))
mypermanova_All_100b <- adonis2(all_Types_100b@otu_table ~ Type, method = "bray", data = data.frame(all_Types_100b@sam_data))# Variety and Biomass

write.xlsx(data.frame(mypermanova_All_100f), file = "C:/Users/User/Desktop/PERMANOVA_Types_Fungi.xlsx")
write.xlsx(data.frame(mypermanova_All_100b), file = "C:/Users/User/Desktop/PERMANOVA_Types_Bacteria.xlsx")

# PERMANOVA for factor = Biomass
mypermanova_All_100f_biomass <- adonis2(all_Types_100f@otu_table ~ Biomass, method = "bray", data = data.frame(all_Types_100f@sam_data))
mypermanova_All_100b_biomass <- adonis2(all_Types_100b@otu_table ~ Biomass, method = "bray", data = data.frame(all_Types_100b@sam_data))# Variety and Biomass

write.xlsx(data.frame(mypermanova_All_100f_biomass), file = "C:/Users/User/Desktop/PERMANOVA_Biomass_Fungi.xlsx")
write.xlsx(data.frame(mypermanova_All_100b_biomass), file = "C:/Users/User/Desktop/PERMANOVA_Biomass_Bacteria.xlsx")

# PERMANOVA for factor = Varieties in each Type
# For fungi
mypermanova_C_100f <- adonis2(fun_C_100@otu_table ~ Variety, method = "bray", data = data.frame(fun_C_100@sam_data))
mypermanova_GU_100f <- adonis2(fun_GU_100@otu_table ~ Variety, method = "bray", data = data.frame(fun_GU_100@sam_data))
mypermanova_GR_100f <- adonis2(fun_GR_100@otu_table ~ Variety, method = "bray", data = data.frame(fun_GR_100@sam_data))

write.xlsx(data.frame(mypermanova_C_100f), file = "C:/Users/User/Desktop/PERMANOVA_Variety_C_Fungi.xlsx")
write.xlsx(data.frame(mypermanova_GU_100f), file = "C:/Users/User/Desktop/PERMANOVA_Variety_GU_Fungi.xlsx")
write.xlsx(data.frame(mypermanova_GR_100f), file = "C:/Users/User/Desktop/PERMANOVA_Variety_GR_Fungi.xlsx")

# For Bacteria
mypermanova_C_100b <- adonis2(bac_C_100@otu_table ~ Variety, method = "bray", data = data.frame(bac_C_100@sam_data))
mypermanova_GU_100b <- adonis2(bac_GU_100@otu_table ~ Variety, method = "bray", data = data.frame(bac_GU_100@sam_data))
mypermanova_GR_100b <- adonis2(bac_GR_100@otu_table ~ Variety, method = "bray", data = data.frame(bac_GR_100@sam_data))

write.xlsx(data.frame(mypermanova_C_100b), file = "C:/Users/User/Desktop/PERMANOVA_Variety_C_Bacteria.xlsx")
write.xlsx(data.frame(mypermanova_GU_100b), file = "C:/Users/User/Desktop/PERMANOVA_Variety_GU_Bacteria.xlsx")
write.xlsx(data.frame(mypermanova_GR_100b), file = "C:/Users/User/Desktop/PERMANOVA_Variety_GR_Bacteria.xlsx")

##-------------------------------------------------------------------------------------------------------
# NMDS Analysis taking account both factors: Types and varieties
# For fungi and yeasts
all_Types_100f <- transform_sample_counts(fungi_rare, function(OTU) 100*OTU/sum(OTU))

ord.nmds.bray1f <- ordinate(all_Types_100f, method="NMDS", distance="bray")

all_Types_100b <- transform_sample_counts(bacteria_rare, function(OTU) 100*OTU/sum(OTU))

ord.nmds.bray1b <- ordinate(all_Types_100b, method="NMDS", distance="bray")

pdf(file = "NMDSv2.pdf",   # The directory you want to save the file in
    width = 10,
    height = 10)

a <- plot_ordination(all_Types_100f, ord.nmds.bray1f, color="Type", label = NULL, title=paste("NMDS Analysis (stress ",round(ord.nmds.bray1f$stress, 2),")", sep = "")) +
  geom_point(aes(shape = Variety), size = 4)+
  scale_shape_manual(values = c(07:25))+
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha=1/4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size = 10))+
  theme(legend.position = "bottom")+
  labs(shape = "Variety\n/Combination")

# For bacteria
b <- plot_ordination(all_Types_100b, ord.nmds.bray1b, color="Type", label = NULL, title=paste("NMDS Analysis (stress ",round(ord.nmds.bray1b$stress, 2),")", sep = "")) +
  geom_point(aes(shape = Variety), size = 4)+
  scale_shape_manual(values = c(07:25))+
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha=1/4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size = 10))+
  theme(legend.position = "bottom")+
  labs(shape = "Variety\n/Combination")

c <- plot_ordination(all_Types_100f, ord.nmds.bray1f, color="Biomass", label = NULL, title=paste("NMDS Analysis (stress ",round(ord.nmds.bray1f$stress, 2),")", sep = "")) +
  #geom_point(aes(shape = Variety), size = 4)+
  scale_shape_manual(values = c(07:25))+
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha=1/4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size = 10))+
  theme(legend.position = "bottom")+
  labs(shape = "Variety\n/Combination")

# For bacteria
d <- plot_ordination(all_Types_100b, ord.nmds.bray1b, color="Biomass", label = NULL, title=paste("NMDS Analysis (stress ",round(ord.nmds.bray1b$stress, 2),")", sep = "")) +
  #geom_point(aes(shape = Variety), size = 4)+
  scale_shape_manual(values = c(07:25))+
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha=1/4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size = 10))+
  theme(legend.position = "bottom")+
  labs(shape = "Variety\n/Combination")

top <- ggarrange(
  a, b,
  ncol = 2,
  nrow = 1,
  align = "hv",
  common.legend = TRUE,
  legend = "bottom",
  labels = c("(a)", "(b)"),
  label.x = 0.08,
  label.y = 0.07,
  font.label = list(
    size = 16,
    face = "bold",
    color = "black"
  )
)

bottom <- ggarrange(
  c, d,
  ncol = 2,
  nrow = 1,
  align = "hv",
  common.legend = TRUE,
  legend = "bottom",
  labels = c("(c)", "(d)"),
  label.x = 0.08,
  label.y = 0.07,
  font.label = list(
    size = 16,
    face = "bold",
    color = "black"
  )
)

final_fig <- ggarrange(
  top,
  bottom,
  ncol = 1,
  nrow = 2,
  align = "v"
)

print(final_fig)

dev.off()

sessionInfo()
