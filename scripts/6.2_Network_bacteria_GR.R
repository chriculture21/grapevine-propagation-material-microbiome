# ============================================================================
# Bacterial co-occurrence network, grafted rooted vines (GR)
# SPRING associations at genus level, nodes coloured by phylum
# Mirrors the GU script so the two figures are directly comparable
# ============================================================================

library(NetCoMi)
library(phyloseq)

out_dir <- "C:/Users/User/Desktop"

tax_b <- data.frame(tax_table(bacteria_cleaned))
xyl <- rownames(tax_b)[grepl("Xylella", tax_b$Genus)]
xyl

tax_table(bacteria_cleaned)[xyl, "Genus"] <- "Xanthomonadaceae_unclassified"
tax_table(bacteria_cleaned)[xyl, ]
# ============================================================================
# 1. Network construction
# ============================================================================
bac_GR    <- subset_samples(bacteria_cleaned, Type == "GR")
gr_genus  <- tax_glom(bac_GR, taxrank = "Genus")

net_spring_GR <- netConstruct(gr_genus,
                              taxRank     = "Genus",
                              filtTax     = "highestFreq",
                              filtTaxPar  = list(highestFreq = 50),
                              filtSamp    = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              measure     = "spring",
                              measurePar  = list(nlambda = 40,
                                                 rep.num = 50,
                                                 Rmethod = "approx"),
                              normMethod  = "none",
                              zeroMethod  = "none",
                              sparsMethod = "none",
                              dissFunc    = "signed",
                              verbose     = 2,
                              seed        = 123456)

props_spring_GR <- netAnalyze(net_spring_GR,
                              centrLCC    = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar      = "eigenvector",
                              weightDeg   = FALSE,
                              normDeg     = FALSE)


# ============================================================================
# 2. Map nodes to phylum
#    Palette is keyed to the GU phylum list so colours match across figures.
# ============================================================================
tax_gr <- data.frame(tax_table(gr_genus))
tax_gr$Genus  <- sub("^[a-z]__", "", tax_gr$Genus)
tax_gr$Phylum <- sub("^[a-z]__", "", tax_gr$Phylum)

nodes_GR       <- colnames(net_spring_GR$adjaMat1)
node_phylum_GR <- tax_gr$Phylum[match(nodes_GR, tax_gr$Genus)]
node_phylum_GR[is.na(node_phylum_GR)] <- "Unclassified"

# fixed phylum-colour key shared by both networks
phy_key <- c("Actinobacteriota"  = "#E69F00",
             "Bacteroidota"      = "#56B4E9",
             "Campilobacterota"  = "#009E73",
             "Desulfobacterota"  = "#CC79A7",
             "Firmicutes"        = "#0072B2",
             "Proteobacteria"    = "#D55E00",
             "Spirochaetota"     = "#F0E442",
             "Verrucomicrobiota" = "#999999")

# any phylum in GR that was not in GU gets a grey slot - check this prints empty
setdiff(unique(node_phylum_GR), names(phy_key))

feat_GR        <- factor(node_phylum_GR, levels = names(phy_key))
feat_GR        <- droplevels(feat_GR)
names(feat_GR) <- nodes_GR

phyla_GR <- levels(feat_GR)
pal_GR   <- unname(phy_key[phyla_GR])

table(node_phylum_GR)


# ============================================================================
# 3. Figure
# ============================================================================
draw_network_GR <- function() {
  plot(props_spring_GR,
       nodeColor      = "feature",
       featVecCol     = feat_GR,
       colorVec       = pal_GR,
       nodeSize       = "eigenvector",
       cexNodes       = 2,
       nodeTransp     = 25,
       hubBorderCol   = "black",
       repulsion      = 0.92,
       rmSingles      = "all",
       labelScale     = FALSE,
       labelFont      = 3,
       hubLabelFont   = 4,
       cexLabels      = 1.5,
       title1         = "Bacterial network on GR",
       showTitle      = TRUE,
       cexTitle       = 2.3,
       mar            = c(2, 3, 4, 14))
  
  legend("topright", inset = c(0.12, 0.1), bty = "n", cex = 2.5,
         title = "Phylum", legend = phyla_GR,
         pch = 21, pt.cex = 2.2, col = "grey30", pt.bg = pal_GR)
  
  legend("topright", inset = c(0.12, 0.60), bty = "n", cex = 2.5,
         title = "Estimated association:",
         legend = c("+", "-"), lty = 1, lwd = 3,
         col = c("#009900", "red"), horiz = TRUE)
  
  usr <- par("usr")
  text(usr[1] + 0.02 * diff(usr[1:2]), usr[3] + 0.03 * diff(usr[3:4]),
       "(b)", font = 2, cex = 1.8)      # (b) - GU figure is (a)
}

pdf(file.path(out_dir, "Network_GR.pdf"), width = 24, height = 14,
    useDingbats = FALSE)
draw_network_GR()
dev.off()

# ============================================================================
# 4. Numbers for the figure legend
# ============================================================================
cat("Samples        :", ncol(net_spring_GR$countMat1), "\n")
cat("Nodes          :", length(nodes_GR), "\n")
cat("Edges          :", sum(net_spring_GR$adjaMat1 != 0) / 2, "\n")
cat("Positive edges :", sum(net_spring_GR$assoMat1 > 0 & net_spring_GR$adjaMat1 != 0) / 2, "\n")
cat("Negative edges :", sum(net_spring_GR$assoMat1 < 0 & net_spring_GR$adjaMat1 != 0) / 2, "\n")

props_spring_GR$globalProps
props_spring_GR$hubs

clust_tab_GR <- data.frame(Genus   = names(props_spring_GR$clustering$clust1),
                           Cluster = as.integer(props_spring_GR$clustering$clust1))
clust_tab_GR$Phylum <- as.character(feat_GR[clust_tab_GR$Genus])
clust_tab_GR <- clust_tab_GR[order(clust_tab_GR$Cluster, clust_tab_GR$Genus), ]

write.csv(clust_tab_GR, file.path(out_dir, "Network_GR_clusters.csv"),
          row.names = FALSE)


# ============================================================================
# 5. Formal GU vs GR comparison
#    Both networks must be built in one netConstruct call on the same taxa.
# ============================================================================
net_comp <- netConstruct(data  = amgut_genus,      # GU (from the GU script)
                         data2 = gr_genus,         # GR
                         taxRank     = "Genus",
                         filtTax     = "highestFreq",
                         filtTaxPar  = list(highestFreq = 50),
                         filtSamp    = "totalReads",
                         filtSampPar = list(totalReads = 1000),
                         measure     = "spring",
                         measurePar  = list(nlambda = 40,
                                            rep.num = 50,
                                            Rmethod = "approx"),
                         normMethod  = "none",
                         zeroMethod  = "none",
                         sparsMethod = "none",
                         dissFunc    = "signed",
                         verbose     = 2,
                         seed        = 123456)

props_comp <- netAnalyze(net_comp,
                         centrLCC    = FALSE,      # FALSE when comparing
                         clustMethod = "cluster_fast_greedy",
                         hubPar      = "eigenvector",
                         weightDeg   = FALSE,
                         normDeg     = FALSE)

# permutation test - slow. Start with 100 permutations to check it runs,
# then rerun with 1000 for the manuscript.
comp <- netCompare(props_comp,
                   permTest = TRUE,
                   nPerm    = 100,
                   cores    = 1,            # Windows: no forking
                   seed     = 123456)

summary(comp, groupNames = c("GU", "GR"))

# side-by-side plot on a shared layout
pdf(file.path(out_dir, "Network_GU_vs_GR.pdf"), width = 24, height = 13)
plot(props_comp,
     sameLayout    = TRUE,
     layoutGroup   = "union",
     rmSingles     = "inboth",
     nodeSize      = "eigenvector",
     labelScale    = FALSE,
     labelFont     = 3,
     cexLabels     = 1.0,
     groupNames    = c("GU - grafted unrooted", "GR - grafted rooted"),
     cexTitle      = 2.2)
dev.off()


# ============================================================================
# Notes
# ----------------------------------------------------------------------------
# * Key question: does the anaerobic module seen in GU persist, fragment, or
#   disappear after six months in the field?
# * netCompare reports differences in global properties (modularity, edge
#   density, average path length) and in centrality per taxon, with Jaccard
#   indices for hub-set overlap and adjusted p-values.
# * Report n for each group - they differ (61 GR vs ~56 GU after filtering).
# ============================================================================
