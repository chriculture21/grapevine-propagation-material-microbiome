library(phyloseq)
library(vegan)

set.seed(1234)

setwd("~path")
bacteria_cleaned <- readRDS("bacteria_cleaned.RDS")

meta_b <- data.frame(sample_data(bacteria_cleaned))

# counts matrix, samples as rows (vegan orientation)
counts_b <- as(otu_table(bacteria_cleaned), "matrix")
if (taxa_are_rows(bacteria_cleaned)) counts_b <- t(counts_b)

lib_b <- rowSums(counts_b)
obs_b <- specnumber(counts_b)


# ============================================================================
# Rarefaction analysis
# ============================================================================
summary(lib_b)
sort(lib_b)[1:10]        # long lower tail; smallest library is 97 reads

rarecurve(counts_b, step = 50, col = "darkblue", label = FALSE, xlim = c(0, 5000),
          xlab = "Reads", ylab = "Observed ASVs",
          main = "Rarefaction curves - 16S")
abline(v = 400, col = "red", lty = 2)
mtext("(b)", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)

# Candidate depths: samples kept, reads used, richness recovered
do.call(rbind, lapply(c(400, 600, 800, 1000, 2000, 3000), function(d) {
  keep <- lib_b >= d
  data.frame(depth            = d,
             samples_kept     = sum(keep),
             pct_reads_used   = round(100 * sum(pmin(lib_b[keep], d)) / sum(lib_b), 1),
             pct_richness_rec = round(100 * mean(rarefy(counts_b[keep, ], d) / obs_b[keep]), 1))
}))


# ============================================================================
# Is low depth confounded with the design?
# ============================================================================
for (d in c(400, 1000, 2000)) {
  cat("\n--- depth", d, "---\n")
  print(table(meta_b$Type, lib_b >= d))
}

kruskal.test(lib_b ~ meta_b$Type)
tapply(lib_b, meta_b$Type, median)

# ============================================================================
# Rarefy
# ============================================================================
bacteria_rare <- rarefy_even_depth(bacteria_cleaned,
                                   sample.size = 400,
                                   rngseed     = 1234,
                                   replace     = FALSE,
                                   trimOTUs    = TRUE)

bacteria_rare

# Which samples were excluded
setdiff(sample_names(bacteria_cleaned), sample_names(bacteria_rare))

saveRDS(bacteria_rare, "~path/bacteria_rarefied.rds")

sessionInfo()
