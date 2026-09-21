library(phyloseq)
library(dplyr)
library(ggplot2)
library(paletteer)
library(ggtext)
library(writexl)

# set working directory
wd <- setwd("~path")

# Load cleaned phyloseq object
ps <- readRDS("fungi_cleaned.rds")
ps1 <- readRDS("bacteria_Cleaned.rds")

# if you want to study a subset of taxa or samples, you have to subset_taxa
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps1 <- prune_samples(sample_sums(ps1) > 0, ps1)

# Convert to relative abundance (%)
ps_rel <- transform_sample_counts(ps, function(x) 100 * x / sum(x))
ps_rel1 <- transform_sample_counts(ps1, function(x) 100 * x / sum(x))

otu_mat <- as(otu_table(ps_rel), "matrix")
otu_mat1 <- as(otu_table(ps_rel1), "matrix")

# Aggregate at Genus level
ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
ps_genus1 <- tax_glom(pse_rel1, taxrank = "Genus", NArm = FALSE)

# Melt to long format
ps_melt <- psmelt(ps_genus)
ps_melt1 <- psmelt(ps_genus1)

# Make sure Genus is character
ps_melt$Genus <- as.character(ps_melt$Genus)
ps_melt1$Genus <- as.character(ps_melt1$Genus)

# Group low-abundance phyla as "< 1%"
ps_melt <- ps_melt %>%
  group_by(Type, Genus) %>%
  mutate(median_abund = median(Abundance)) %>%
  ungroup()

keep_taxa <- unique(ps_melt$Genus[ps_melt$median_abund > 1])
ps_melt$Genus[!(ps_melt$Genus %in% keep_taxa)] <- "< 1%"

ps_melt1 <- ps_melt1 %>%
  group_by(Type, Genus) %>%
  mutate(median_abund = median(Abundance)) %>%
  ungroup()

keep_taxa1 <- unique(ps_melt1$Genus[ps_melt1$median_abund > 2.5])
ps_melt1$Genus[!(ps_melt1$Genus %in% keep_taxa1)] <- "< 2.5%"

# Summarise mean abundance per Type and Phylum
ps_melt_sum <- ps_melt %>%
  group_by(Type, Variety, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(Type, Variety) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance)) %>%  # rescale to 100%
  ungroup()

ps_melt_sum$Genus <- sub("^[a-z]__", "", ps_melt_sum$Genus) # for fungi

ps_melt_sum$Type <- factor(
  ps_melt_sum$Type,
  levels = c("C", "GU", "GR")
)

ps_melt_sum1 <- ps_melt1 %>%
  group_by(Type, Variety, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(Type, Variety) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance)) %>%  # rescale to 100%
  ungroup()

ps_melt_sum1$Type <- factor(
  ps_melt_sum1$Type,
  levels = c("C", "GU", "GR")
)

italicize_legend <- function(x) {
  
  higher <- "(mycota|mycotina|mycetes|mycetidae|ales|ineae|aceae|oideae|eae)$"
  
  lapply(x, function(tax) {
    
    # Long Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium name
    if (grepl(
      "Allorhizobium.*Neorhizobium.*Pararhizobium.*Rhizobium",
      tax
    )) {
      
      # Split into two lines and italicize both
      bquote(
        atop(
          italic("Allorhizobium−Neorhizobium"),
          italic("-Pararhizobium−Rhizobium")
        )
      )
      
    } else if (
      tax == "< 1%" ||
      grepl(higher, tax, ignore.case = TRUE) ||
      grepl(
        "unidentified|unclassified|incertae",
        tax,
        ignore.case = TRUE
      )
    ) {
      
      # Higher taxa and <1% remain roman
      bquote(.(tax))
      
    } else {
      
      # Genus names italicized
      bquote(italic(.(tax)))
    }
  })
}

# Plot stacked bar chart
pdf(file = "Barplotsv2.pdf",   # The directory you want to save the file in
    width = 18,
    height = 22)

a <- ggplot(ps_melt_sum, aes(x = Type, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.8, color="black") +
  facet_wrap(~ Type, scales = "free_x")+
  labs(x = NULL, y = NULL, title = "Fungi and yeasts", fill = "Taxon") +
  theme_classic(base_size = 20) +
  scale_fill_paletteer_d(
    "ggthemes::calc",
    labels = italicize_legend
  )+
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),,
    legend.title = element_text(size = 22),
    legend.box.margin = margin(l = -15),
    legend.key.height = grid::unit(0.7, "cm"),
    legend.key.width = grid::unit(0.7, "cm"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5))

b <- ggplot(ps_melt_sum1, aes(x = Type, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  facet_wrap(~ Type, scales = "free_x")+
  labs(x = NULL, y = NULL, title = "Bacteria", fill = "Taxon") +
  theme_classic(base_size = 20) +
  scale_fill_paletteer_d(
    "ggthemes::calc",
    labels = italicize_legend
  )+
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.box.margin = margin(l = -15),
    legend.key.height = grid::unit(0.7, "cm"),
    legend.key.width = grid::unit(0.7, "cm"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5))

fig <- ggarrange(a, b, ncol = 1, nrow = 2,
                 align = "h", common.legend = FALSE,
                 labels = c("(a)", "(b)"),
                 label.x = 0.02, label.y = 0.05,
                 font.label = list(size = 30, face = "bold", color = "black"))

print(
  annotate_figure(fig,
                  bottom = text_grob("Variety/Rootstock", size = 30, face = "bold"),
                  left   = text_grob("Relative abundance (%)", size = 30, face = "bold", rot = 90))
)
  

dev.off()

# save the files with the abundances
write_xlsx(ps_melt_sum, 'fungi_abundances.xlsx')
write_xlsx(ps_melt_sum1, 'bacteria_abundances.xlsx')

sessionInfo()
