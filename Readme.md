# 🧬 miRNA Differential Expression & Enrichment App

This interactive Shiny application allows users to perform and explore **differential expression analysis of miRNA-seq data** using DESeq2, visualize results through various dimensionality reduction and plotting tools, perform **functional enrichment**, run **Random Forest classification**, and conduct **power analysis** to support experimental design.

---

## ✨ Features

- 📈 DESeq2-based differential expression analysis
- 🌋 Volcano plot, PCA, UMAP, heatmaps, and barplots
- 🧠 Functional enrichment (Enrichr & clusterProfiler fallback)
- 🌲 Random Forest classification with AUC, sensitivity, and feature importance
- 🔋 Power analysis for sample size estimation across miRNAs
- 💾 Exportable tables and plots
- 📖 Built-in README tab for user guidance
- 🧊 Dark theme UI with Skyrim-style styling (optional)

---

## 🛠️ Tech Stack

| Layer              | Tool/Library            |
|-------------------|-------------------------|
| Language           | R (≥ 4.0.0)             |
| Frontend           | Shiny, Plotly, DT       |
| DE Analysis        | DESeq2                  |
| Enrichment         | enrichR, clusterProfiler, org.Hs.eg.db |
| ML Classification  | randomForest, pROC      |
| Power Analysis     | pwr                     |
| Containerization   | Docker, Singularity/Apptainer |
| Optional Styling   | Custom CSS              |

---

## 🚀 Installation

### 📦 Prerequisites

- [R (≥ 4.0)](https://www.r-project.org/)
- R packages: `shiny`, `DESeq2`, `ggplot2`, `dplyr`, `plotly`, `umap`, `pheatmap`, `enrichR`, `clusterProfiler`, `org.Hs.eg.db`, `randomForest`, `pROC`, `pwr`, `DT`
- [Docker](https://docs.docker.com/) or [Apptainer/Singularity](https://apptainer.org/)
- Optional: RStudio or R command line for local use

---

## 🧪 Local Use in R

```r
# Clone repository and run from R or RStudio
shiny::runApp("path/to/app")

Ensure the following files exist:

app.R

/data/mirna_targets.csv

/www/styles.css

readme.txt

🐳 Docker Deployment
Build
in bash
docker build -t mirna-de-app .

Then open: http://localhost:3838 in your browser.


🧬 Singularity / Apptainer Deployment (for HPC)

Step 1: Build the container
in bash
apptainer build mirna-de-app.sif Singularity.def
Step 2: Run interactively on an HPC login node (for testing)
in bash
apptainer shell mirna-de-app.sif
cd /srv/shiny-server
shiny-server
Step 3: Submit to a compute node (e.g. via SLURM)
Create a job script run_app.sbatch:

in bash
#!/bin/bash
#SBATCH --job-name=mirna_app
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=shiny_app.out

module load apptainer
apptainer run mirna-de-app.sif
Then submit with


sbatch run_app.sbatch
Check the output log for the URL to connect (you may need to expose ports through your cluster's web gateway or use port forwarding with SSH).


📂 Data Structure
project/
├── app.R
├── Singularity.def
├── Dockerfile
├── data/
│   └── mirna_targets.csv
├── www/
│   └── styles.css
├── readme.txt

📖 Example Input Requirements
Counts File: A CSV matrix of raw counts (rows = miRNAs, columns = samples)

Metadata File: CSV file with sample names and conditions (last column should be the condition)



Built with ❤️ and R to make bioinformatics just a little more magical.





