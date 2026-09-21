library(phyloseq)
library(vegan)

set.seed(1234)

setwd("~path")
fungi_cleaned <- readRDS("fungi_cleaned.RDS")

meta <- data.frame(sample_data(fungi_cleaned))

# Counts matrix, samples as rows (vegan orientation)
counts <- as(otu_table(fungi_cleaned), "matrix")
if (taxa_are_rows(fungi_cleaned)) counts <- t(counts)

lib <- rowSums(counts)
obs <- specnumber(counts)


# ============================================================================
# Rarefaction analysis
# ============================================================================
summary(lib)
sort(lib)[1:10]

rarecurve(counts, step = 100, col = "darkblue", label = FALSE, xlim = c(0, 12000),
          xlab = "Reads", ylab = "Observed OTUs",
          main = "Rarefaction curves - ITS2")
abline(v = min(lib), col = "red", lty = 2)
mtext("(a)", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)

# Candidate depths: samples kept, reads used, richness recovered
do.call(rbind, lapply(c(1206, 2000, 3000, 4000, 5000, 8000, 10000), function(d) {
  keep <- lib >= d
  data.frame(depth            = d,
             samples_kept     = sum(keep),
             pct_reads_used   = round(100 * sum(pmin(lib[keep], d)) / sum(lib), 1),
             pct_richness_rec = round(100 * mean(rarefy(counts[keep, ], d) / obs[keep]), 1))
}))


# ============================================================================
# Is low depth confounded with the design?
# ============================================================================
table(meta$Type,    lib >= 3000)      # all shallow libraries are GU
table(meta$Biomass, lib >= 3000)      # balanced

kruskal.test(lib ~ meta$Type)
tapply(lib, meta$Type, median)


# ============================================================================
# Rarefy
# ============================================================================
fungi_rare <- rarefy_even_depth(fungi_cleaned,
                                sample.size = min(lib),
                                rngseed     = 1236,
                                replace     = FALSE,
                                trimOTUs    = TRUE)

fungi_rare        # 1530 taxa, 154 samples; 276 OTUs dropped

saveRDS(fungi_rare, "C:/Users/User/Desktop/fungi_rarefied.rds")

sessionInfo()

