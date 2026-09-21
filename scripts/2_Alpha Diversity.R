## Alpha diversity indices for fungi and bacteria between the different types of propagation material
##C: Canes, GU: Grafted unrooted vines, GR: Grafted rooted vines

library(vegan)
library(phyloseq)
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(openxlsx)

# set working directory
setwd("~path")

# load rds files after rarefraction analysis
fungi_rare <- readRDS("fungi_rarefied.rds")
bacteria_rare <- readRDS("bacteria_rarefied.rds")

fungi_richness <- phyloseq::estimate_richness(fungi_rare, measures = c("Observed", "Shannon"))
fungi_richness <- cbind(fungi_richness, data.frame(sample_data(fungi_rare)))
bacteria_richness <- phyloseq::estimate_richness(bacteria_rare, measures = c("Observed", "Shannon"))
bacteria_richness <- cbind(bacteria_richness, data.frame(sample_data(bacteria_rare)))

setwd('~path')
write.xlsx(fungi_richness, file = "Fungi_richness.xlsx", rowNames = TRUE)
write.xlsx(bacteria_richness, file = "Bacteria_richness.xlsx", rowNames = TRUE)

head(fungi_richness)
head(bacteria_richness)

# create a list of pairwise comaprisons
ps1.meta <- meta(fungi_rare)
bmi <- levels(ps1.meta$Type) # get the variables

# make a pairwise list for comparison
pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])

print(pairs)

# Plot alpha diversity richness estimates
pdf(file = "C:/Users/User/Desktop/Alpha diversity.pdf",
    width = 12,
    height = 15)

a <- plot_richness(fungi_rare, x='Type', measures = c("Observed","Shannon"))+
  geom_violin(aes(fill=Type), stat = 'ydensity')+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  scale_x_discrete(limits = c("C","GU","GR")
  )+
  theme(text = element_text(size=22))+
  theme(legend.position = "bottom")+
  stat_compare_means(comparisons = pairs)+
  labs(
    x = NULL,
    y = NULL,
    title = "Fungi and yeasts")

b <- plot_richness(bacteria_rare, x='Type', measures = c("Observed","Shannon"))+
  geom_violin(aes(fill=Type), stat = 'ydensity')+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  scale_x_discrete(limits = c("C","GU","GR")
  )+
  theme(text = element_text(size=22))+
  theme(legend.position = "bottom")+
  stat_compare_means(comparisons = pairs)+
  labs(
    x = NULL,
    y = NULL,
    title = "Bacteria")

fig <- ggarrange(a, b, ncol = 1, nrow = 2,
          align = "v", common.legend = TRUE,
          labels = c("(a)", "(b)"),
          label.x = 0.02, label.y = 0.05,
          font.label = list(size = 30, face = "bold", color = "black"))

print(
  annotate_figure(fig,
                  bottom = text_grob("Propagation material type", size = 30, face = "bold"),
                  left   = text_grob("Alpha diversity index", size = 30, face = "bold", rot = 90))
)

dev.off()

sessionInfo()
