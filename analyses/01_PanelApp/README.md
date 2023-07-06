# 01. Kidney-Genetics Analysis: PanelApp

This repository contains an R script that retrieves and processes genetic panel data related to kidney diseases from the PanelApp project in the UK and Australia. The script is part of a broader project called "kidney-genetics".

## Requirements

- R (version 3.6.0 or above)
- R libraries: readr, tidyverse, rvest, jsonlite, config

## Usage

Clone the repository and run the R script located in `analyses/01_PanelApp/`. 

```bash
git clone https://github.com/halbritter-lab/kidney-genetics
cd kidney-genetics/analyses/01_PanelApp/
Rscript 01_PanelApp.R
```

Before running the script, set the `CONFIG_FILE` environment variable to the path of your config file.

```bash
export CONFIG_FILE="/path/to/your/config/file"
```

## Description

The script performs the following actions:

1. Sets up the working directory based on the configuration file.
2. Loads HGNC (HUGO Gene Nomenclature Committee) functions from an external R script.
3. Retrieves kidney disease related panels from both PanelApp UK and PanelApp Australia, meaning all panels that include "renal" or "kidney" in its name.
4. Combines all kidney related PanelApp panels into a single tibble (data frame).
5. Processes the panels to create a unified format for all genes.
6. Saves the processed data into CSV files. 

The output of the script is two CSV files:

- `01_PanelApp_panels.<date>.csv`: Contains the data for each kidney-related panel retrieved from PanelApp.
- `01_PanelApp_genes.<date>.csv`: Contains the data for each gene related to kidney diseases, processed and formatted for further analysis.

## License

This project is licensed under the terms of the MIT license. See the LICENSE file for details.
