# 02. Kidney-Genetics Analysis: Literature

This repository contains an R script for analyzing kidney genetics literature. The script downloads and processes data from various publications, and generates a summary of the analyzed genes. We assesed the publications using following search query: 
        (1)	"Kidney"[Mesh] OR "Kidney Diseases"[Mesh] OR kidney OR renal
	AND (2)	"Genetic Structures"[Mesh] OR "Genes"[Mesh] OR genetic test OR gene panel OR gene panels OR multigene panel OR multigene panels OR targeted panel*

## Prerequisites

Make sure you have the following libraries installed:

- `readxl`: Used to load Excel files.
- `tidyverse`: Required for data processing.
- `jsonlite`: Needed for HGNC functions.
- `officer`: Used for docx files.
- `pdftools`: Required for pdf files.
- `config`: Needed for config loading.

For this analysis you need to have access (VPN) to to respective journals to download the supplement files.

## Usage

Clone the repository and run the R script located in `analyses/02_Literature/`. 

```bash
git clone https://github.com/halbritter-lab/kidney-genetics
cd kidney-genetics/analyses/02_Literature/
Rscript 02_Literature.R
```

## Configuration

Before running the script, set the `CONFIG_FILE` environment variable to the path of your config file.

```bash
export CONFIG_FILE="/path/to/your/config/file"
```

## Description

The script performs the following actions:

1. Load the required libraries.
2. Load global functions from the `hgnc-functions.R` file.
3. Load the manually curated overview Excel table and filter the data.
4. Download all files specified in the Excel table.
5. Load and normalize different files from publications.
6. Generate a summary of the analyzed genes.
7. The results will be saved in the `results` directory with a timestamp appended to the filename.

## License

This project is licensed under the [MIT License](LICENSE).