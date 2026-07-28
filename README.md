# Inferring Predation and Locomotor Strategies

## Overview

This repository supports the project **Convergent Carnivores or Divergent Dasyurids? Inferring Predation and Locomotor Strategies of Extant and Extinct Carnivorous Marsupials from Locomotor Shape**.

The study examines whether locomotor shape is associated with predatory strategy across carnivoran and dasyuromorph mammals. The analyses combine side-view stride sequences, two-dimensional landmark digitisation, geometric morphometrics, principal component analysis, ecological clustering, phylogenetically informed analyses, trajectory analysis, and classification of the thylacine.

The final dataset contains:

- **49 species**
- **115 included stride sequences**
- **74 symmetrical strides**
- **41 asymmetrical strides**

Only sequences marked **Checked** in `Morpho_Meta_V4_URLs.xlsx` are included in the analyses and in Table S1.

---

## Repository structure

| Folder / file | Description |
|---|---|
| `data/` | Raw and processed stride coordinates and analysis-ready datasets. |
| `figures/` | Figures generated for the manuscript and supplementary material. |
| `Morpho_Meta_V4_URLs.xlsx` | Video-level metadata, including inclusion status, source URL or provenance, video identifier, gait, symmetry class, species, mass estimate, and screening comments. |
| `Table_S1_Analysed_Strides.docx` | Supplementary Table S1. Lists all included stride sequences with species information, gait counts, video identifiers, URLs or provenance statements, and comments. |
| `Centroid_analysis_Symgaits.R` | Centroid and morphospace analyses for symmetrical gaits. |
| `Centroid_analysis_aSymgaits.R` | Centroid and morphospace analyses for asymmetrical gaits. |
| `Correlations_angles_along_PCs.R` | Correlates anatomical angles and limb ratios with principal-component axes. |
| `TrajectoryAnalysis.R` | Compares trajectory size, shape, and orientation among hunting clusters. |
| `thylacine_cluster_classification.R.R` | Classifies the thylacine centroid among hunting-strategy clusters. |
| `thylacine_cluster_classification_figure.R` | Produces the thylacine classification figure. |
| `Make_Supp_Tables.R` | Generates supplementary statistical tables. |
| `Make_density_plots_and_Stickfigures.R` | Produces density plots and stick-figure visualisations. |
| `PLotting_wireframes_alongPCs.R` | Produces wireframes illustrating shape variation along principal-component axes. |
| `README.md` | Repository documentation. |

---

## Video metadata and provenance

`Morpho_Meta_V4_URLs.xlsx` contains the video-level metadata used to document source provenance and screening decisions. Public URLs are supplied where available. Other records are identified as author-owned footage or by a source/provenance note when a direct public URL was not recorded.

The comments fields preserve notes made during video screening, including departures from an ideal perpendicular view, turning near the end of a stride, head position, discontinuities, individual identities in release footage, and source-attribution notes. Table S1 reports only the 124 rows marked **Checked**; excluded candidate clips remain in the metadata workbook for transparency but are not used in the analyses.

### Table S1 - species, gait, and source information

Table S1 reports, for every included stride:

- scientific and common name;
- clade and hunting-strategy cluster;
- species-level counts of symmetrical and asymmetrical strides;
- video identifier and gait;
- source URL or provenance statement; and
- screening comments.

For the complete list, see [`Table_S1_Analysed_Strides.docx`](Table_S1_Analysed_Strides.docx).

---

## Getting started

### Prerequisites

The analyses were developed in R. The required packages include:

```r
install.packages(c(
  "geomorph", "vegan", "lme4", "car", "effectsize", "performance",
  "RRPP", "phytools", "dendextend", "ggplot2", "dplyr", "tidyr",
  "cowplot", "patchwork"
))
```

### Clone the repository

```bash
git clone https://github.com/cclemente/Inferring-Predation-and-Locomotor-Strategies-.git
cd Inferring-Predation-and-Locomotor-Strategies-
```

### Run the analyses

The scripts address complementary parts of the analysis rather than forming a single automated pipeline. File paths may need to be updated for the local clone.

1. Run the symmetrical and asymmetrical centroid analyses.
2. Run the angle-to-PC correlation analyses.
3. Run the trajectory analyses.
4. Run the thylacine classification and figure scripts.
5. Generate the density plots, wireframes, figures, and supplementary tables.

---

## Methods summary

- **Video selection:** full strides from healthy adult animals moving approximately perpendicular to the camera on level ground.
- **Digitisation:** 19 two-dimensional landmarks were digitised across each stride.
- **Temporal standardisation:** each stride was interpolated to 100 frames.
- **Shape analysis:** landmark configurations were Procrustes-aligned and analysed using principal component analysis.
- **Ecological clustering:** species were grouped using six dimensions of hunting ecology.
- **Statistical analyses:** PERMANOVA, linear mixed-effects models, phylogenetically informed analyses, and trajectory comparisons were used to test locomotor differences.
- **Thylacine inference:** the thylacine's walking centroid was compared with hunting-strategy cluster centroids using multiple classifiers.

---

## Data and source-video use

The repository provides derived coordinate data, metadata, scripts, and figures used for reproducibility. Rights to externally hosted source videos remain with their respective creators or owners. URLs and provenance statements are provided to document the material used in the analysis.

---

## Citation

Please cite the associated manuscript when it becomes available. Until then, the repository may be cited as:

> Gaschk, J. L., Cieri, R. L., Chaseling, B. R., Hamilton, D. G., & Clemente, C. J. (2026). *Inferring Predation and Locomotor Strategies* [GitHub repository]. https://github.com/cclemente/Inferring-Predation-and-Locomotor-Strategies-

---

## Authors

Joshua L. Gaschk  
Robert L. Cieri  
Bella R. Chaseling  
David G. Hamilton  
Christofer J. Clemente

[Biomechanics and Biorobotics Lab](https://biomechanics-and-biorobotics-lab.com/), University of the Sunshine Coast, Australia


Last updated: 27 July 2026
