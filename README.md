# 🐾 Inferring Predation and Locomotor Strategies

## Overview  
This repository supports the project exploring **how locomotor behavior reflects predatory strategy** across carnivoran and dasyurid species.  
The analysis combines high-speed video sequences, geometric morphometrics, and gait classification to infer functional and ecological groupings from movement.

We extracted and analyzed **124 side-view stride sequences** from **84 individuals (40 asymmetric gaits)** spanning multiple clades and behavioral clusters (e.g. *Social Hunters*, *Opp. Pouncers*, *Opp. Grapplers*, *Large Grapplers*).  
Key outputs include centroid trajectory analyses, wireframe visualizations of limb movement, and supplementary tables for publication.

---

## 🗂 Repository Structure

| Folder / File | Description |
|----------------|-------------|
| `data/` | Contains raw and processed stride data, shape coordinates, and metadata. |
| `Centroid_analysis_Symgaits.R` | R script for centroid trajectory analysis of symmetrical gaits. |
| `Centroid_analysis_aSymgaits.R` | R script for centroid trajectory analysis of asymmetrical gaits. |
| `Correlations_angles_along_PCs.R` | Correlates limb/track angles with principal component axes. |
| `LAD_attempt2.R`, `LAD_attempt2_figure.R` | Implements Linear Additive Decomposition (LAD) analysis and plots. |
| `Make_Supp_Tables.R` | Generates all supplementary tables (e.g. Table S1). |
| `Make_density_plots_and_Stickfigures.R` | Produces density plots of shape variation and stick-figure visualisations. |
| `PLotting_wireframes_alongPCs.R` | Plots wireframes of limb/track shapes along PC axes. |
| `TrajectoryAnalysis.R` | Performs stride trajectory analysis and kinematic summary metrics. |
| `Table_S1_Species_List.docx` | Supplementary species list and video metadata (Table S1). |
| `README.md` | This document. |

---

## 📊 Data Description

### Table S1 – Species and Gait Summary  
This table lists each species included in the analysis, their common names, taxonomic clade, behavioral cluster, and counts of symmetrical versus asymmetrical gait sequences.  
It also provides the specific video files used for each animal.  

- **Total species:** 84  
- **Symmetrical gaits:** 84  
- **Asymmetrical gaits:** 40  
- **Total sequences:** 124  

Example entry:

| Species | Common name | Clade | Cluster | Symm. | Asymm. | Videos |
|:---------|:-------------|:------|:---------|------:|------:|:--------|
| *Acinonyx jubatus* | Cheetah | Felid | Social Hunter | 2 | 4 | Cheetah2, Cheetah3, CheetahDog1, Cheetah1 |

For the full list, see [`Table_S1_Species_List.docx`](./Table_S1_Species_List.docx).

---

## 🧰 Getting Started

### Prerequisites  
You’ll need **R (≥ 4.0)** and several common R packages:  
```r
install.packages(c("geomorph", "ggplot2", "dplyr", "tidyr", "cowplot", "patchwork"))
```

### Clone the Repository  
```bash
git clone https://github.com/cclemente/Inferring-Predation-and-Locomotor-Strategies.git
cd Inferring-Predation-and-Locomotor-Strategies
```

### Run the Analyses  
Open R or RStudio and run scripts in sequence:  
1. `Centroid_analysis_Symgaits.R` → analyze symmetrical gait trajectories  
2. `Centroid_analysis_aSymgaits.R` → analyze asymmetrical gait trajectories  
3. `Correlations_angles_along_PCs.R` → correlate shape features with PCA axes  
4. `Make_density_plots_and_Stickfigures.R` → create density and stick-figure plots  
5. `Make_Supp_Tables.R` → regenerate supplementary tables  

Figures and tables will be saved to the current working directory or to folders specified in each script.

---

## 🧬 Methods Summary  
- **Data acquisition:** high-speed video footage from public datasets and field recordings.  
- **Digitization:** stride outlines were digitized and Procrustes-aligned using `geomorph`.  
- **PCA:** performed on aligned shapes to extract major modes of stride variation.  
- **Behavioral clustering:** species classified as *Social Hunters*, *Opp. Pouncers*, *Opp. Grapplers*, *Large Grapplers*, *Anteaters*, etc.  
- **Analysis outputs:** centroid trajectories, gait symmetry indices, and visual morphospaces.  

---

## 🤝 Contributing  
Contributions are welcome!  
If you’d like to add species, stride sequences, or additional analyses:

1. Fork this repository  
2. Add your files to the appropriate `data/` subfolder  
3. Update metadata (e.g. `species_metadata.csv`)  
4. Re-run `Make_Supp_Tables.R` to integrate the new entries  
5. Submit a pull request with a concise description of the update

---

## 📜 Citation  
If you use this dataset or scripts, please cite:

> Clemente, C. J. (2025). *Inferring Predation and Locomotor Strategies.*  
> GitHub Repository: [https://github.com/cclemente/Inferring-Predation-and-Locomotor-Strategies](https://github.com/cclemente/Inferring-Predation-and-Locomotor-Strategies)

---

## ⚖️ License  
This project is distributed under the **MIT License**.  
See the [`LICENSE`](./LICENSE) file for details.

---

## 🧠 Author  
**Christofer J. Clemente**  
Biomechanics & Biorobotics Lab,  
University of the Sunshine Coast, Australia  
[https://www.christoferclemente.com](https://www.christoferclemente.com)

---

*Last updated: 23 October 2025*
