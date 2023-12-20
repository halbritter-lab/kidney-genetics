# B. Annotation HGNC

This script (`B_AnnotationHGNC.R`) is designed to annotate HGNC data within the nephrology-focused kidney-genetics project. The script fetches data from various databases like HGNC, OMIM, gnomAD, clinvar, and genCC. It then processes and saves the results for subsequent analyses.

## Prerequisites

Before running the script, make sure you have the following libraries installed:

- `tidyverse`  (for general table operations)
- `biomaRt`    (to get gene coordinates)
- `STRINGdb`   (to compute StringDB identifiers)
- `R.utils`    (to gzip downloaded and result files)
- `readxl`     (to read xlsx files)
- `config`     (for config loading)
- `janitor`    (for cleaning column names)

## Usage

Run the R script located in the directory where it resides. You may need to adjust paths and configurations according to your specific setup.

```bash
Rscript B_AnnotationHGNC.R
```

Before running the script, set the `CONFIG_FILE` environment variable to the path of your config file.

```bash
export CONFIG_FILE="/path/to/your/config/file"
```

## Script Overview

The script performs the following steps:

1. Loads the required libraries.
2. Defines project-related variables (topic, name, script path).
3. Reads configuration variables and sets the working directory.
4. Sets specific global options.
5. Loads global functions from external scripts.
6. Downloads required database sources from HGNC, OMIM, gnomAD, clinvar, and genCC.

The results are saved and ready for further analysis or visualization.

### License

This project is licensed under the terms of the MIT license. See the LICENSE file for details.
