# Fungal and bacterial communities across grapevine propagation material

R scripts for the downstream analysis of ITS2 (fungi and yeasts) and 16S rRNA (bacteria) amplicon sequencing data from different types of grapevine propagation material. The scripts cover rarefaction, alpha and beta diversity, taxonomic composition, differential abundance, and bacterial co-occurrence networks.

> **Associated publication:** Manuscript in preparation

## Study design

Samples were grouped by three factors, which correspond to columns in the phyloseq `sample_data`:

| Variable | Levels | Description |
|---|---|---|
| `Type` | `C`, `GU`, `GR` | Type of propagation material: C = Canes, GU = grafted unrooted vines, GR = grafted rooted vines|
| `Variety` | Variety, rootstock, or scion/rootstock combination |
| `Biomass` `Low`, `High` | Low, or high pathogen biomass based on qPCR assays |

## Workflow

Raw reads were processed with DADA2 (quality control, filtering, denoising, and taxonomic assignment) before these scripts were run. The scripts start from the two cleaned phyloseq objects, `fungi_cleaned.RDS` and `bacteria_Cleaned.RDS`.

Alpha and beta diversity were calculated on rarefied data. Taxonomic composition, differential abundance, and network analyses used the un-rarefied cleaned data.

```
DADA2 pre-processing (not included)
        │
        ▼
fungi_cleaned.RDS / bacteria_Cleaned.RDS
        │
        ├── 1.2_Rarefaction_analysis_fungi.R ──┐
        ├── 1.1_Rarefaction_analysis_bacteria.R ┤──► *_rarefied.rds
        │                                   │         │
        │                                   │         ├── 2_Alpha_Diversity.R
        │                                   │         └── 3_Beta Diversity.R
        │
        ├── 4_Composition.R
        ├── 5_ANCOMBC2_DA.R
        ├── 6.1_Network_bacteria_GU.R
        └── 6.2_Network_bacteria_GR.R
```

## Scripts

### 1. Rarefaction

**`1.1_Rarefaction_analysis_fungi.R`** and **`1.2_Rarefaction_analysis_bacteria.R`**

These scripts plot rarefaction curves, tabulate candidate rarefaction depths (samples retained, reads used, and richness recovered at each depth), and test whether sequencing depth is confounded with `Type` using a Kruskal–Wallis test. The rarefaction depth for each dataset was chosen to retain as much of the richness as possible without removing samples unevenly across the experimental design.

Rarefaction used `phyloseq::rarefy_even_depth()` without replacement. Outputs are `fungi_rarefied.rds` and `bacteria_rarefied.rds`.

Singletons had already been removed during cleaning, so singleton-based richness estimators (e.g. Chao1) are not used.

### 2. Alpha diversity

**`2_Alpha_Diversity.R`**

Calculates Observed richness and Shannon diversity on the rarefied data, compares Types pairwise with Wilcoxon tests (`ggpubr::stat_compare_means`), and plots violin/box plots for fungi and bacteria.

Outputs: `Fungi_richness.xlsx`, `Bacteria_richness.xlsx`, `Alpha diversity.pdf`

### 3. Beta diversity

**`3_Beta Diversity.R`**

Uses Bray–Curtis dissimilarities on relative abundances calculated from the rarefied data.

- PERMANOVA (`vegan::adonis2`) testing the effects of `Type` and `Biomass` across all samples, and of `Variety` within each Type
- Pairwise PERMANOVA (`pairwiseAdonis::pairwise.adonis`, Benjamini–Hochberg correction) between Types, and between Varieties within each Type
- NMDS ordinations coloured by `Type` and by `Biomass`

Outputs: one `.xlsx` table per test, `NMDSv2.pdf`

### 4. Taxonomic composition

**`4_Composition.R`**

Converts counts to relative abundance, agglomerates to genus level, and plots stacked bar charts of mean relative abundance per Type and Variety. Genera below a median relative abundance of 1% (fungi) or 2.5% (bacteria) are grouped together.

Outputs: `Barplotsv2.pdf`, `fungi_abundances.xlsx`, `bacteria_abundances.xlsx`

### 5. Differential abundance

**`5_ANCOMBC2_DA.R`**

Runs ANCOM-BC2 at genus level on the unrarefied data, with `Type` as the fixed effect and C as the reference level, so the contrasts are C vs GU and C vs GR. Settings: minimum library size 1,000 reads, prevalence filter 10%, FDR correction.

Genera were retained at q < 0.05 for fungi and q < 0.001 for bacteria, and classified as enriched or depleted using |log fold change| ≥ 1.

Outputs: `DA_Fungi.pdf`, `DA_bacteria.pdf`, and one `.xlsx` table per contrast

### 6. Bacterial co-occurrence networks

**`6.1_Network_bacteria_GU.R`** and **`6.2_Network_bacteria_GR.R`**

Genus-level co-occurrence networks for grafted unrooted (GU) and grafted rooted (GR) material, built with NetCoMi using SPRING associations on the 50 most frequent genera (StARS model selection: 40 λ values, 50 subsamples; samples with fewer than 1,000 reads excluded). Nodes are sized by eigenvector centrality and coloured by phylum, and hubs are identified by eigenvector centrality. Clusters were detected with the fast-greedy algorithm.

- The GU script also tests whether network clusters are structured by oxygen requirement (Fisher's exact test). This requires adding an `Oxygen` column (aerobe / anaerobe / facultative, assigned from the literature) to the exported cluster table.
- The GR script uses the same phylum colour key as the GU script so the figures are directly comparable. It then builds both networks together and compares them with `netCompare()` (permutation test).
- In the GR script, ASVs assigned to *Xylella* are relabelled as `Xanthomonadaceae_unclassified` [state the reason, e.g. assignment not supported by ...].

Outputs: `Network_GU.pdf`, `Network_GR.pdf`, `Network_GU_vs_GR.pdf`, `Network_GU_clusters.csv`, `Network_GR_clusters.csv`

## Requirements

R 4.4.3 with the following packages:

| Source | Packages |
|---|---|
| CRAN | vegan, dplyr, tibble, tidyr, ggplot2, forcats, stringr, ggpubr, ggtext, paletteer, openxlsx, writexl, knitr |
| Bioconductor | phyloseq, microbiome, ANCOMBC |
| GitHub | [pairwiseAdonis](https://github.com/pmartinezarbizu/pairwiseAdonis), [NetCoMi](https://github.com/stefpeschel/NetCoMi) |

## Usage

1. Place `fungi_cleaned.RDS` and `bacteria_Cleaned.RDS` in a `data/` folder. These files are not yet included in the repository.
2. Update the file paths at the top of each script if needed.
3. Run the rarefaction scripts first, then the remaining scripts in the order shown in the workflow above. Run `Network_bacteria_GU.R` before `Network_bacteria_GR.R`, because the network comparison uses objects created in the GU script.

Seeds are set for rarefaction and network construction.

## Data availability

The associated manuscript is in preparation. The cleaned phyloseq objects (`fungi_cleaned.RDS`, `bacteria_Cleaned.RDS`) and the accession numbers for the raw sequencing reads will be released here after publication. Until then, the scripts are shared to document the analysis workflow.

## License

[MIT / GPL-3.0 / CC BY 4.0]

## Contact

Christos D. Tsoukas, [chriculture@gmail.com / ORCID: https://orcid.org/0009-0008-0628-9542]
