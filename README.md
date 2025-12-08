# DIAToolComp

# DIATool — Unified DIA Proteomics Data Converter

Welcome to **DIATool**! This tool simplifies the comparison and integration of results from different **Data-Independent Acquisition (DIA)** proteomics software by converting their output into a **standardized format**.

We support conversion from popular DIA analysis platforms including:
- **EncyclopeDIA**
- **DIA-NN**
- **Skyline**
- **OpenSWATH**
- **Spectronaut**

Whether you're benchmarking tools, performing meta-analyses, or building downstream pipelines, DIATool ensures your data speaks the same language. 🧬📊

🔗 **Try it live**: [https://www.guomics.com/DIAtool/](https://www.guomics.com/DIAtool/)

---

## 📁 Repository Overview

- **`app.R`**  
  The core **Shiny web application** that powers the interactive interface. It handles file uploads, format detection, data conversion, and result downloads.

- **`template.RData`**  
  A sample dataset used internally by the app for demonstration and testing purposes. It includes example matrices in the unified format to illustrate expected outputs.

- **`ReportConversion/`**  
  Contains backend conversion scripts and logic for parsing and standardizing outputs from each supported DIA tool.

- **`ResultPlot/`**  
  Utilities for visualizing converted results (e.g., QC metrics, overlap plots)—used in the online version for exploratory analysis.

---

## 🚀 How It Works

1. Upload your result file(s) from any supported DIA software.
2. DIATool automatically detects the source tool and converts the data into a consistent tabular format at three levels:
    - **Precursor**
    - **Peptide**
    - **Protein**
3. Download the standardized matrices for downstream analysis or cross-tool comparison.

This unified schema enables seamless integration with statistical tools, visualization libraries, or machine learning workflows.

---

## 💻 Run Locally (Optional)

To run the app on your own machine:

```r
# Install required packages (if needed)
install.packages(c("shiny", "dplyr", "readr", "stringr"))

# Clone the repo
git clone https://github.com/guomics-lab/DIAToolComp.git
cd DIAToolComp

# Launch the app
Rscript -e "shiny::runApp()"