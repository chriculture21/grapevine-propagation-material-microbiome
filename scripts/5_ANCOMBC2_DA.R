# Differential abundance analysis between the different types of propagation material with ANCOMBC2
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
ps <- readRDS("fungi_cleaned.rds") #not rarefied
ps2 <- readRDS("bacteria_cleaned.rds") #not rarefied

# Run differential abundance for fungi
out <- ancombc2(
  data = ps,            # phyloseq object
  tax_level = "Genus",
  fix_formula = "Type + Biomass",
  struc_zero = TRUE,
  neg_lb = FALSE,
  verbose = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE,
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  trend_control = NULL,
  group = "Type",             # your main variable
  rand_formula = NULL,
  alpha = 0.05,
  p_adj_method = "holm",
  lib_cut = 1206,             # min library size
  prv_cut = 0.10,             # prevalence filter (10%)
  s0_perc = 0.05,
  n_cl = 2,                   # number of threads
  global = TRUE
)

# and for bacteria
out2 <- ancombc2(
  data = ps2,
  tax_level = "Genus",
  fix_formula = "Type + Biomass",
  struc_zero = TRUE,
  neg_lb = FALSE,
  verbose = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE,
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  trend_control = NULL,
  group = "Type",             # main variable
  rand_formula = NULL,
  alpha = 0.05,
  p_adj_method = "holm",
  lib_cut = 400,             # min library size
  prv_cut = 0.10,             # prevalence filter (10%)
  s0_perc = 0.05,
  n_cl = 2,                   # number of threads
  global = TRUE
)

# Get the results table for both fungi and bacteria
results <- out$res #for fungi
results2 <- out2$res #for bacteria

# remove prefix from fungal taxonomy
results <- results %>%
  mutate(taxon = gsub("^[a-z]__", "", taxon))

# Select only the desired columns from the results and generate the data frames
sig_ancombc1 <- results %>%
  filter(q_TypeGU < 0.05, abs(lfc_TypeGU) > 1) %>%
  dplyr::select(taxon, lfc_TypeGU, se_TypeGU, q_TypeGU)

sig_ancombc2 <- results %>%
  filter(q_TypeGR < 0.05, abs(lfc_TypeGR) > 1) %>%
  dplyr::select(taxon, lfc_TypeGR, se_TypeGR, q_TypeGR)

sig_ancombc3 <- results2 %>%
  filter(q_TypeGU < 0.001, abs(lfc_TypeGU) > 1) %>%
  dplyr::select(taxon, lfc_TypeGU, se_TypeGU, q_TypeGU)

sig_ancombc4 <- results2 %>%
  filter(q_TypeGR < 0.001, abs(lfc_TypeGR) > 1) %>%
  dplyr::select(taxon, lfc_TypeGR, se_TypeGR, q_TypeGR)

# Add a new column in the dataframe, indicating if a genus is enriched or depleted q-value and log fold change
# C vs GU - Fungi
sig_ancombc1 <- sig_ancombc1 %>%
  mutate(Significance = case_when(
    lfc_TypeGU >=1 & q_TypeGU <= 0.05 ~ 'Enriched',
    lfc_TypeGU <=-1 & q_TypeGU <= 0.05 ~ 'Depleted',
    abs(lfc_TypeGU) < 1 | q_TypeGU > 0.05 ~ 'Not significant'
  ))

sig_ancombc1 = sig_ancombc1 %>%
  arrange(lfc_TypeGU)

# C vs GR - Fungi
sig_ancombc2 <- sig_ancombc2 %>%
  mutate(Significance = case_when(
    lfc_TypeGR >=1 & q_TypeGR <= 0.05 ~ 'Enriched',
    lfc_TypeGR <=-1 & q_TypeGR <= 0.05 ~ 'Depleted',
    abs(lfc_TypeGR) < 1 | q_TypeGR > 0.05 ~ 'Not significant'
  ))

sig_ancombc2 = sig_ancombc2 %>%
  arrange(lfc_TypeGR)

# C vs GU - Bacteria
sig_ancombc3 <- sig_ancombc3 %>%
  mutate(Significance = case_when(
    lfc_TypeGU >=1 & q_TypeGU <= 0.001 ~ 'Enriched',
    lfc_TypeGU <=-1 & q_TypeGU <= 0.001 ~ 'Depleted',
    abs(lfc_TypeGU) < 1 | q_TypeGU > 0.001 ~ 'Not significant',
  ))

sig_ancombc3 = sig_ancombc3 %>%
  arrange(lfc_TypeGU)

# C vs GR - Bacteria
sig_ancombc4 <- sig_ancombc4 %>%
  mutate(Significance = case_when(
    lfc_TypeGR >=1 & q_TypeGR <= 0.05 ~ 'Enriched',
    lfc_TypeGR <=-1 & q_TypeGR <= 0.05 ~ 'Depleted',
    abs(lfc_TypeGR) < 1 | q_TypeGR > 0.05 ~ 'Not significant'
  ))

sig_ancombc4 = sig_ancombc4 %>%
  arrange(lfc_TypeGR)

# Only Genus names in italics for the plot
italicize_genus <- function(x) {
  higher <- "(mycota|mycotina|mycetes|mycetidae|ales|ineae|aceae|oideae|eae)$"
  vapply(x, function(lab) {
    p <- strsplit(lab, " \\| ")[[1]]
    prefix <- if (length(p) >= 2) paste0(p[1], " | ") else ""
    tax <- p[length(p)]                                  # last part = taxon name
    if (grepl(higher, tax) ||
        grepl("unidentified|unclassified|incertae", tax, ignore.case = TRUE)) {
      lab                                                # leave roman
    } else {
      paste0(prefix, "*", tax, "*")                      # genus → italic
    }
  }, character(1), USE.NAMES = FALSE)
}

## Create lfc barplots with standard error - Fungi
pdf(file = "~DA_Fungi.pdf",   # The directory you want to save the file in
    width = 14,
    height = 8)

a <- ggplot(sig_ancombc1, aes(x = taxon, y = lfc_TypeGU, fill = Significance)) +
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
    subtitle = "q-value cutoff = 0.05; lfc cutoff = 1"
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

b <- ggplot(sig_ancombc2, aes(x = taxon, y = lfc_TypeGR, fill = Significance)) +
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
    subtitle = "q-value cutoff = 0.05; lfc cutoff = 1"
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
                 legend = "top",
                 labels = c("(a)", "(b)"),
                 label.x = 0, label.y = 1,
                 font.label = list(size = 20, face = "bold", color = "black")) #arrange the two plots side by side

print(
  annotate_figure(fig,
                  bottom = text_grob("Log Fold Change (LFC)", size = 18, face = "bold"),
                  left   = text_grob("Taxon", size = 18, face = "bold", rot = 90))
)

dev.off() #close the devise

# lfc barplots for bacteria
n_tax  <- max(nrow(sig_ancombc3), nrow(sig_ancombc4))
fs     <- 8

italicize_genus <- function(x) {
  higher <- "(mycota|mycotina|mycetes|mycetidae|ales|ineae|aceae|oideae|eae)$"
  out <- vapply(x, function(lab) {
    p <- strsplit(lab, " \\| ")[[1]]
    prefix <- if (length(p) >= 2) paste0(p[1], " | ") else ""
    tax <- p[length(p)]
    if (grepl(higher, tax) ||
        grepl("unidentified|unclassified|incertae", tax, ignore.case = TRUE)) {
      lab
    } else {
      paste0(prefix, "*", tax, "*")
    }
  }, character(1), USE.NAMES = FALSE)
  
  # wrap the long Rhizobium-group name over two lines
  sub("Allorhizobium-Neorhizobium-", "Allorhizobium-Neorhizobium-*<br>*", out, fixed = TRUE)
}

pdf(file = "~DA_bacteria.pdf",   # The directory you want to save the file in
    width = 14,
    height = 8)

c <- ggplot(sig_ancombc3, aes(reorder(taxon, lfc_TypeGU), y = lfc_TypeGU, fill = Significance)) +
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
    subtitle = "q-value cutoff = 0.001; lfc cutoff = 1"
  ) +
  theme_gray() +
  theme(
    axis.text.y = element_markdown(size = 14, face = 'italic'),
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
    legend.position = "none"
  )

d <- ggplot(sig_ancombc4, aes(x = reorder(taxon, lfc_TypeGR), y = lfc_TypeGR, fill = Significance)) +
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
    subtitle = "q-value cutoff = 0.001; lfc cutoff = 1"
  ) +
  theme_gray() +
  theme(
    axis.text.y = element_markdown(size = 14, face = 'italic'),
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
    legend.position = "none"
  )

fig <- ggarrange(c, d, ncol = 2, nrow = 1,
                 align = "hv", common.legend = FALSE, legend = NULL,
                 labels = c("(c)", "(d)"),
                 label.x = 0, label.y = 1,
                 font.label = list(size = 20, face = "bold", color = "black")) #arrange the two plots side by side

print(
  annotate_figure(fig,
                  bottom = text_grob("Log Fold Change (LFC)", size = 18, face = "bold"),
                  left   = text_grob("Taxon", size = 18, face = "bold", rot = 90))
)

dev.off() # close device

# Save the data in xlsx files
write.xlsx(sig_ancombc1, "~DA_CvsGU_Fungi.xlsx")
write.xlsx(sig_ancombc2, "~DA_CvsGR_Fungi.xlsx")
write.xlsx(results, "~All_DA_Fungi.xlsx")
write.xlsx(sig_ancombc3, "~DA_CvsGU_Bacteria.xlsx")
write.xlsx(sig_ancombc4, "~DA_CvsGR_Bacteria.xlsx")
write.xlsx(results2, "~All_DA_Bacteria.xlsx")

# Save session
sessionInfo()
