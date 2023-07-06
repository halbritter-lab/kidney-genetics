# 03. Kidney-Genetics Analysis: Diagnostic Panels

This repository contains a script that performs data extraction and analysis for Kidney-Genetics Diagnostic Panels. The script downloads data from different web URLs, identifies and extracts genes related to kidney diseases, standardizes the gene list per panel, and then summarizes the data. Following Diagnostic Panels are included:

1. [Centogene Nephrology](https://www.centogene.com/diagnostics/our-tests/ngs-panels/nephrology)
2. [CeGaT kidney diseases](https://cegat.com/diagnostics/rare-diseases/kidney-diseases/)
3. [Preventiongenetics comprehensive inherited kidney diseases](https://www.preventiongenetics.com/testInfo?val=Comprehensive-Inherited-Kidney-Diseases-Panel)
4. [Invitae progressive renal diseases](https://www.invitae.com/en/providers/test-catalog/test-75000)
5. [Invitae expanded renal diseases](https://www.invitae.com/en/providers/test-catalog/test-633100)
6. [MGZ Nephrology](https://www.mgz-muenchen.de/gen-panels/panel/nierenerkrankungen-und-elektrolyte-gesamtpanel.html)
7. [MVZ kidney diseases](https://www.medizinische-genetik.de/diagnostik/humangenetik/panelinformation?panel=229)
8. [Natera renasight comprehensive kidney genes](https://www.natera.com/organ-health/renasight-genetic-testing/gene-conditions-list/)
9. [Mayocliniclabs renal genetics](https://www.mayocliniclabs.com/test-catalog/Overview/618086)
10. [Blueprintgenetics Nephrology](https://blueprintgenetics.com/tests/panels/nephrology/)

## Required Libraries

This script uses the following R libraries:

- `readr` and `readxl` for reading data.
- `tidyverse` for data manipulation and visualization.
- `rvest` and `jsonlite` for web scraping and handling JSON data.
- `curl` and `httr` for handling URLs and HTTP.
- `config` for loading configuration variables.
- `webdriver` for headless browsing.

Please install these packages in R before running the script.
You will also require PhantomJS installed on your system.

## Usage

Clone the repository and run the R script located in `analyses/03_DiagnosticPanels/`. 

```bash
git clone https://github.com/halbritter-lab/kidney-genetics
cd kidney-genetics/analyses/03_DiagnosticPanels/
Rscript 03_DiagnosticPanels.R
```

## Configuration

Before running the script, set the `CONFIG_FILE` environment variable to the path of your config file.

```bash
export CONFIG_FILE="/path/to/your/config/file"
```

## Description

The script performs the following actions:

1. Setting up the environment and loading necessary libraries.
2. Reading in the list of diagnostic panels to be used.
3. Downloading the relevant data using PhantomJS.
4. Extracting the list of genes associated with each diagnostic panel.
5. Summarizing and writing out the results.


## Attention

The current functionality is focused on scraping data from specific web sources and requires manual intervention when these sources change or update their website structure. The script may not perform as expected if the web pages' structure is significantly different from when the script was written.
For example the `xpath` argument used in several of the `html_nodes` calls specifies the location of the data in the web page's HTML structure. This might need to be modified if the structure of the web pages changes.

## License

The scripts and documentation in this project are released under the [MIT License](LICENSE).