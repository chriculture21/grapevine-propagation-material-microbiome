# ============================================================================
# Bacterial co-occurrence network, grafted unrooted vines (GU)
# SPRING associations at genus level, nodes coloured by phylum
# ============================================================================

library(NetCoMi)
library(phyloseq)

setwd("~path")
bacteria_cleanes <- readRDS("bacteria_Cleaned.rds")

# Network construction
bac_GU      <- subset_samples(bacteria_cleaned, Type == "GU")
amgut_genus <- tax_glom(bac_GU, taxrank = "Genus")

net_spring <- netConstruct(amgut_genus,
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

props_spring <- netAnalyze(net_spring,
                           centrLCC    = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar      = "eigenvector",
                           weightDeg   = FALSE,
                           normDeg     = FALSE)


# Map nodes to phylum
tax <- data.frame(tax_table(amgut_genus))
tax$Genus  <- sub("^[a-z]__", "", tax$Genus)
tax$Phylum <- sub("^[a-z]__", "", tax$Phylum)

nodes       <- colnames(net_spring$adjaMat1)
node_phylum <- tax$Phylum[match(nodes, tax$Genus)]
node_phylum[is.na(node_phylum)] <- "Unclassified"

feat        <- factor(node_phylum)
names(feat) <- nodes                       # names are essential here

phyla <- levels(feat)
pal   <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7",
           "#0072B2", "#D55E00", "#F0E442", "#999999")[seq_along(phyla)]

table(node_phylum)
stopifnot(length(pal) == length(phyla))    # extend pal if you have >8 phyla


# Figure
draw_network <- function() {
  plot(props_spring,
       nodeColor      = "feature",         # NOT "cluster" - featVecCol is
       featVecCol     = feat,              # matched by node name
       colorVec       = pal,               # one colour per phylum level
       nodeSize       = "eigenvector",
       cexNodes       = 2,
       nodeTransp     = 25,
       hubBorderCol   = "black",
       repulsion      = 0.92,              # spreads nodes, reduces overlap
       rmSingles      = "all",             # drop unconnected genera
       labelScale     = FALSE,
       labelFont      = 3,                 # italic genus names
       hubLabelFont   = 4,                 # bold italic hubs
       cexLabels      = 1.5,
       title1         = "Bacterial network on GU",
       showTitle      = TRUE,
       cexTitle       = 2.3,
       mar            = c(2, 3, 4, 14))    # right margin holds the legends
  
  legend("topright", inset = c(0.12, 0.1), bty = "n", cex = 2.5,
         title = "Phylum", legend = phyla,
         pch = 21, pt.cex = 2.2, col = "grey30", pt.bg = pal)
  
  legend("topright", inset = c(0.12, 0.60), bty = "n", cex = 2.5,
         title = "Estimated association:",
         legend = c("+", "-"), lty = 1, lwd = 3,
         col = c("#009900", "red"), horiz = TRUE)
  
  usr <- par("usr")
  text(usr[1] + 0.02 * diff(usr[1:2]), usr[3] + 0.03 * diff(usr[3:4]),
       "(a)", font = 2, cex = 1.8)
}

pdf(file.path(out_dir, "Network_GU.pdf"), width = 24, height = 14,
    useDingbats = FALSE)
draw_network()
dev.off()


# Numbers for the figure legend and methods
cat("Samples        :", ncol(net_spring$countMat1), "\n")
cat("Nodes          :", length(nodes), "\n")
cat("Edges          :", sum(net_spring$adjaMat1 != 0) / 2, "\n")
cat("Positive edges :", sum(net_spring$assoMat1 > 0 & net_spring$adjaMat1 != 0) / 2, "\n")
cat("Negative edges :", sum(net_spring$assoMat1 < 0 & net_spring$adjaMat1 != 0) / 2, "\n")

props_spring$globalProps
props_spring$hubs



# Cluster membership (the anaerobic-guild question)
clust_tab <- data.frame(Genus   = names(props_spring$clustering$clust1),
                        Cluster = as.integer(props_spring$clustering$clust1))
clust_tab$Phylum <- as.character(feat[clust_tab$Genus])
clust_tab <- clust_tab[order(clust_tab$Cluster, clust_tab$Genus), ]

write.csv(clust_tab, file.path(out_dir, "Network_GU_clusters.csv"),
          row.names = FALSE)

clust_tab <- read.xlsx("Network_GU_clusters.xlsx")

# Add an Oxygen column (aerobe / anaerobe / facultative) from the literature,
# then test whether clusters are structured by oxygen requirement:
fisher.test(table(clust_tab$Cluster, clust_tab$Oxygen), simulate.p.value = TRUE, B = 10000 )

table(clust_tab$Cluster, clust_tab$Oxygen)

sessionInfo()