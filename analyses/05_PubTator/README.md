# 05. Kidney-Genetics Analysis: PubTator

This repository contains an R script for performing an analysis on kidney genetics data using the PubTator API and other relevant functions. The script retrieves information on genes and their associated publications related to kidney diseases or renal diseases.

## Prerequisites

Before running the script, make sure you have the following libraries installed:

- `readr`
- `tidyverse`
- `httr`
- `jsonlite`
- `config`

Please install these packages in R before running the script.

## Usage

Clone the repository and run the R script located in `analyses/05_PubTator/`. 

```bash
git clone https://github.com/halbritter-lab/kidney-genetics
cd kidney-genetics/analyses/05_PubTator/
Rscript 05_PubTator.R
```

## Script Overview

The script performs the following steps:

1. Loads the required libraries.
2. Sets the working directory and global options.
3. Retrieves the NCBI API key from the configuration file.
4. Loads the necessary functions from external source files.
5. Defines the search query for kidney diseases and relevant keywords.
6. Retrieves the number of pages for the search query.
7. Iterates over the pages and retrieves the gene information using the PubTator API.
8. Filters and processes the retrieved data to extract relevant information.
9. Normalizes the data using specific functions.
10. Saves the results as a CSV file.

For this analysis you need your personal NCBI API key.

The script saves the results as a CSV file in the `results` directory. The file name includes the creation date to differentiate between different analysis runs.

## License

This project is licensed under the [MIT License](LICENSE).