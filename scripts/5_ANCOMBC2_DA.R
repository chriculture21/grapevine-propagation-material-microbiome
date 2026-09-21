# Differential abundance analysis between the different types of propagation material
# with ANCOMBC2
# Load libraries
library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(ggpubr)
library(dplyr)
library(ggtext)
library(openxlsx)

setwd('~path')

# Load data
fungi_cleaned <- readRDS("fungi_cleaned.rds")
bacteria_cleaned <- readRDS("bacteria_Cleaned.rds")

# Rename the phyloseq object for clarity
ps <- fungi_cleaned
ps1 <- bacteria_cleaned

# Aggregate to the desired level (for instance Genus)
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus2 <- tax_glom(ps1, taxrank = "Genus")

# Run differential abundance analysis using ANCOMBC2 for fungi
out <- ancombc2(
  data = ps_genus,            # phyloseq object
  tax_level = NULL,           # NULL if already glommed
  fix_formula = "Type",
  group = "Type",             # your main variable
  rand_formula = NULL,        # random effects if needed
  alpha = 0.05,
  p_adj_method = "fdr",
  lib_cut = 1000,             # min library size
  prv_cut = 0.10,             # prevalence filter (10%)
  s0_perc = 0.05,
  n_cl = 1,                   # number of threads
  global = TRUE
)

# and for bacteria
out2 <- ancombc2(
  data = ps_genus2,           # phyloseq object directly
  tax_level = NULL,           # NULL if already glommed
  fix_formula = "Type",
  group = "Type",             # main variable
  rand_formula = NULL,        # random effects if needed
  alpha = 0.05,
  p_adj_method = "fdr",
  lib_cut = 1000,             # min library size
  prv_cut = 0.10,             # prevalence filter (10%)
  s0_perc = 0.05,
  n_cl = 1,                   # number of threads
  global = TRUE
)

# Get the results table for both fungi and bacteria
results <- out$res
results2 <- out2$res

# Extract the taxonomy from the phyloseq object
tax_df <- as.data.frame(tax_table(ps_genus)) # for fungi
tax_df$taxon <- rownames(tax_df)

tax_df2 <- as.data.frame(tax_table(ps_genus2)) # for bacteria
tax_df2$taxon <- rownames(tax_df2)

# join ANCOMBC results with taxonomy
res_tax <- results %>%
  left_join(tax_df %>% dplyr::select(taxon, Genus), by = "taxon") %>%
  mutate(Genus = gsub("^[a-z]__", "", Genus))

res_tax2 <- results2 %>%
  left_join(tax_df2 %>% dplyr::select(taxon, Genus), by = "taxon") %>%
  mutate(Genus = gsub("^[a-z]__", "", Genus))

# Filter significant genera (adjusted p < 0.05 or < 0.001 in any contrast)
sig_ancombc1 <- res_tax %>%
  filter(q_TypeGU < 0.05) %>%
  dplyr::select(taxon, Genus, lfc_TypeGU, se_TypeGU, p_TypeGU)

sig_ancombc2 <- res_tax %>%
  filter(q_TypeGR < 0.05) %>%
  dplyr::select(taxon, Genus, lfc_TypeGR, se_TypeGR, p_TypeGR)

sig_ancombc3 <- res_tax2 %>%
  filter(q_TypeGU < 0.001) %>%
  dplyr::select(taxon, Genus, lfc_TypeGU, se_TypeGU, p_TypeGU)

sig_ancombc4 <- res_tax2 %>%
  filter(q_TypeGR < 0.001) %>%
  dplyr::select(taxon, Genus, lfc_TypeGR, se_TypeGR, p_TypeGR)

# Add a new column in the dataframe, indicating if a genus is enriched, 
# depleted or not significantly different based on the desired p-value and log fold change
# C vs GU - Fungi
sig_ancombc1 <- sig_ancombc1 %>%
  mutate(Significance = case_when(
    lfc_TypeGU >=1 & p_TypeGU <= 0.05 ~ 'Enriched',
    lfc_TypeGU <=-1 & p_TypeGU <= 0.05 ~ 'Depleted',
    abs(lfc_TypeGU) < 1 | p_TypeGU > 0.05 ~ 'Not Significant'
  ))

sig_ancombc1 = sig_ancombc1 %>%
  mutate(y_label = paste0(taxon, " | ", Genus)) %>%
  arrange(lfc_TypeGU)

sig_ancombc1$y_label <- factor(sig_ancombc1$y_label, levels = sig_ancombc1$y_label)

# C vs GR - Fungi
sig_ancombc2 <- sig_ancombc2 %>%
  mutate(Significance = case_when(
    lfc_TypeGR >=1 & p_TypeGR <= 0.05 ~ 'Enriched',
    lfc_TypeGR <=-1 & p_TypeGR <= 0.05 ~ 'Depleted',
    abs(lfc_TypeGR) < 1 | p_TypeGR > 0.05 ~ 'Not significant'
  ))

sig_ancombc2 = sig_ancombc2 %>%
  mutate(y_label = paste0(taxon, " | ", Genus)) %>%
  arrange(lfc_TypeGR)

sig_ancombc2$y_label <- factor(sig_ancombc2$y_label, levels = sig_ancombc2$y_label)

# C vs GU - Bacteria
sig_ancombc3 <- sig_ancombc3 %>%
  mutate(Significance = case_when(
    lfc_TypeGU >=1 & p_TypeGU <= 0.001 ~ 'Enriched',
    lfc_TypeGU <=-1 & p_TypeGU <= 0.001 ~ 'Depleted',
    abs(lfc_TypeGU) < 1 | p_TypeGU > 0.001 ~ 'Not significant',
  ))

sig_ancombc3 = sig_ancombc3 %>%
  mutate(y_label = paste0(taxon, " | ", Genus)) %>%
  arrange(lfc_TypeGU)

sig_ancombc3$y_label <- factor(sig_ancombc3$y_label, levels = sig_ancombc3$y_label)

# C vs GR - Bacteria
sig_ancombc4 <- sig_ancombc4 %>%
  mutate(Significance = case_when(
    lfc_TypeGR >=1 & p_TypeGR <= 0.001 ~ 'Enriched',
    lfc_TypeGR <=-1 & p_TypeGR <= 0.001 ~ 'Depleted',
    abs(lfc_TypeGR) < 1 | p_TypeGR > 0.001 ~ 'Not significant'
  ))

sig_ancombc4 = sig_ancombc4 %>%
  mutate(y_label = paste0(taxon, " | ", Genus)) %>%
  arrange(lfc_TypeGR)

sig_ancombc4$y_label <- factor(sig_ancombc4$y_label, levels = sig_ancombc4$y_label)

# Genus names in italics for the plot
italicize_genus <- function(x) {
  higher <- "(mycota|mycotina|mycetes|mycetidae|ales|ineae|aceae|oideae|eae)$"
  vapply(x, function(lab) {
    p <- strsplit(lab, " \\| ")[[1]]
    if (length(p) < 2) return(lab)                    # no taxon part
    tax <- p[2]
    if (grepl(higher, tax) ||
        grepl("unidentified|unclassified|incertae", tax, ignore.case = TRUE)) {
      lab                                             # leave roman
    } else {
      paste0(p[1], " | *", tax, "*")                  # genus → italic
    }
  }, character(1), USE.NAMES = FALSE)
}

## Create lfc barplots with standard error - Fungi
pdf(file = "DA_Fungi.pdf",   # The directory you want to save the file in
    width = 14,
    height = 8)

a <- ggplot(sig_ancombc1, aes(x = y_label, y = lfc_TypeGU, fill = Significance)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  geom_errorbar(aes(ymin = lfc_TypeGU - se_TypeGU, ymax = lfc_TypeGU + se_TypeGU), width = 0.2, color = "black") +
  geom_hline(yintercept = 0.0, linewidth=0.7, color = "red", linetype = "dashed")+
  geom_hline(yintercept=c(-1, 1), col="red", linetype = 'dashed')+
  scale_fill_manual(values = c('darkred','darkblue','lightgreen'))+
  scale_x_discrete(labels = italicize_genus) +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    title = "C vs GU",
    subtitle = "p-value cutoff = 0.05; lfc cutoff = 1"
  ) +
  theme_gray() +
  theme(
    axis.text.y = element_markdown(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom"
  )

b <- ggplot(sig_ancombc2, aes(x = y_label, y = lfc_TypeGR, fill = Significance)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  geom_errorbar(aes(ymin = lfc_TypeGR - se_TypeGR, ymax = lfc_TypeGR + se_TypeGR), width = 0.2, color = "black") +
  geom_hline(yintercept = 0.0, linewidth=0.7, color = "red", linetype = "dashed")+
  geom_hline(yintercept=c(-1, 1), col="red", linetype = 'dashed')+
  scale_fill_manual(values = c('darkred','darkblue','lightgreen'))+
  scale_x_discrete(labels = italicize_genus) +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    title = "C vs GR",
    subtitle = "p-value cutoff = 0.05; lfc cutoff = 1"
  ) +
  theme_gray() +
  theme(
    axis.text.y = element_markdown(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom"
  )

fig <- ggarrange(a, b, ncol = 2, nrow = 1,
                 align = "hv", common.legend = TRUE,
                 labels = c("(a)", "(b)"),
                 label.x = 0.02, label.y = 0.05,
                 font.label = list(size = 20, face = "bold", color = "black"))

print(
  annotate_figure(fig,
                  bottom = text_grob("Log Fold Change (LFC)", size = 18, face = "bold"),
                  left   = text_grob("Taxon", size = 18, face = "bold", rot = 90))
)

dev.off()

# lfc barplots for bacteria
n_tax  <- max(nrow(sig_ancombc3), nrow(sig_ancombc4))
fs     <- 8

pdf(file = "DA_bacteria.pdf",   # The directory you want to save the file in
    width = 22,
    height = n_tax * fs * 1.3 / 72 + 3)

c <- ggplot(sig_ancombc3, aes(x = y_label, y = lfc_TypeGU, fill = Significance)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  geom_errorbar(aes(ymin = lfc_TypeGU - se_TypeGU, ymax = lfc_TypeGU + se_TypeGU), width = 0.2, color = "black") +
  geom_hline(yintercept = 0.0, linewidth=0.7, color = "red", linetype = "dashed")+
  geom_hline(yintercept=c(-1, 1), col="red", linetype = 'dashed')+
  scale_fill_manual(values = c('darkred','darkblue','lightgreen'))+
  scale_x_discrete(labels = ~ sub("\\| (.*)$", "| *\\1*", .x)) +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    title = "C vs GU",
    subtitle = "p-value cutoff = 0.001; lfc cutoff = 1"
  ) +
  theme_gray() +
  theme(
    axis.text.y = element_markdown(size = 14),
    axis.text.x = element_text(size = 16),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "italic"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom"
  )

d <- ggplot(sig_ancombc4, aes(x = y_label, y = lfc_TypeGR, fill = Significance)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  geom_errorbar(aes(ymin = lfc_TypeGR - se_TypeGR, ymax = lfc_TypeGR + se_TypeGR), width = 0.2, color = "black") +
  geom_hline(yintercept = 0.0, linewidth=0.7, color = "red", linetype = "dashed")+
  geom_hline(yintercept=c(-1, 1), col="red", linetype = 'dashed')+
  scale_fill_manual(values = c('darkred','darkblue','lightgreen'))+
  scale_x_discrete(labels = ~ sub("\\| (.*)$", "| *\\1*", .x)) +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    title = "C vs GR",
    subtitle = "p-value cutoff = 0.001; lfc cutoff = 1"
  ) +
  theme_gray() +
  theme(
    axis.text.y = element_markdown(size = 14),
    axis.text.x = element_text(size = 16),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "italic"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom"
  )

fig <- ggarrange(c, d, ncol = 2, nrow = 1,
                 align = "h", common.legend = TRUE,
                 labels = c("(a)", "(b)"),
                 label.x = 0.02, label.y = 0.05,
                 font.label = list(size = 20, face = "bold", color = "black"))

print(
  annotate_figure(fig,
                  bottom = text_grob("Log Fold Change (LFC)", size = 18, face = "bold"),
                  left   = text_grob("Taxon", size = 18, face = "bold", rot = 90))
)

dev.off() # close device

# Save the data used for plotting
write.xlsx(sig_ancombc1, "C:/Users/User/Desktop/DA_Fungi_CvsGU.xlsx")
write.xlsx(sig_ancombc2, "C:/Users/User/Desktop/DA_Fungi_CvsGR.xlsx")
write.xlsx(sig_ancombc3, "C:/Users/User/Desktop/DA_Bacteria_CvsGU.xlsx")
write.xlsx(sig_ancombc4, "C:/Users/User/Desktop/DA_Bacteria_CvsGR.xlsx")

# Save session
sessionInfo()
