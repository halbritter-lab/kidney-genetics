# Kidney-Genetics: Merge all analyses sources

This script is used to merge all multiple analysis files and generate a wide table with summarized information.

## Prerequisites

Make sure you have the following libraries installed:

- `readr`: Used to load CSV files.
- `tidyverse`: Required for data processing.
- `tools`: Needed for md5 checksum calculation.


## Usage

Clone the repository and run the R script located in `analyses/`. 

```bash
git clone https://github.com/halbritter-lab/kidney-genetics
cd kidney-genetics/
Rscript MergeAnalysesSources.R
```

## Configuration

Before running the script, set the `CONFIG_FILE` environment variable to the path of your config file.

```bash
export CONFIG_FILE="/path/to/your/config/file"
```

## Description

1. **Load Libraries**: The required libraries (`tidyverse`, `readr`, `tools`) are loaded.

2. **Define Relative Script Path**: The `project_name` and `script_path` variables are set to specify the relative path to the script within the project.

3. **Read Config**: The script reads the configuration variables from the provided configuration file.

4. **Set Working Directory**: The working directory is set to the specified project's directory concatenated with the script's path.

5. **Set Global Options**: Global options are set to suppress scientific notation.

6. **Load and Transform Analysis Files**: The script loads and transforms the analysis files located in various folders (`01_PanelApp/results/`, `02_Literature/results/`, `03_DiagnosticPanels/results/`, `04_HPO/results/`, `05_PubTator/results/`). It filters the files to select only those containing "genes" in their paths. The newest file is selected for each analysis. The script then loads the CSV files and combines them into a single table.

7. **Generate Wide Table**: The script generates a wide table by summarizing the information from the loaded files. It computes `evidence_count` as the sum of lists where the `source_evidence` is TRUE and `list_count` as the sum of all lists where the gene is found (`source_evidence` is TRUE or FALSE). The resulting table is stored in the `results_genes_wider` variable.

8. **Save Results**: The merged table is saved as a CSV file in the `merged` directory with a file name based on the current date and time.