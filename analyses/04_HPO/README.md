# 04. Kidney-Genetics Analysis: HPO term detection in clinical databases

This is an R script designed to perform an analysis on the Human Phenotype Ontology (HPO) and Online Mendelian Inheritance in Man (OMIM) database to identify the genes associated with abnormal kidney conditions. This script queries the ontology and database for specific conditions and file formats and processes the data to produce a list of relevant genes and their associated diseases.

### Requirements
The following R libraries are required:

- tidyverse
- jsonlite
- rvest
- httr
- kableExtra
- config

Please install these packages in R using the install.packages() function.

## Usage

Clone the repository and run the R script located in `analyses/04_HPO-detection/`. 

```bash
git clone https://github.com/halbritter-lab/kidney-genetics
cd kidney-genetics/analyses/04_HPO-detection/
Rscript 04_HPO-detection.R
```

Before running the script, set the `CONFIG_FILE` environment variable to the path of your config file.

```bash
export CONFIG_FILE="/path/to/your/config/file"
```

## Description

The script performs the following actions:

1. Fetch HPO terms that are descendants of 'Abnormality of the upper urinary tract' (HP:0010935).
2. Download phenotype data from the HPO and genemap data from OMIM.
3. Merge into a list of genes and related diseases.
4. Write these lists to CSV files in a results directory.

The script generates two output files, stored in the 'results' directory:

1. A list of all HPO terms that are descendants of the term 'Abnormality of the upper urinary tract'.
2. A list of genes associated with the HPO terms, along with their reported gene name, HGNC id, database id, source evidence and source count.

### License

This project is licensed under the terms of the MIT license. See the LICENSE file for details.