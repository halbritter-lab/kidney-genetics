
# C. Annotate Merged Table

This script (`C_AnnotateMergedTable.R`) is part of the nephrology-focused kidney-genetics project and is designed to annotate the merged tables with data needed for manual curation.

## Prerequisites

Before running the script, ensure you have the following libraries installed:

- `tidyverse`          (for general table operations)
- `jsonlite`           (for HGNC requests)
- `rvest`              (for scraping)
- `readr`              (to read files)
- `tools`              (for checksums)
- `R.utils`            (to gzip downloaded and result files)
- `config`             (to read config files)
- `ontologyIndex`      (to read ontology files)

## Usage

Run the R script located in the directory where it resides. Adjust paths and configurations based on your specific setup.

```bash
Rscript C_AnnotateMergedTable.R
```

Before executing the script, set the `CONFIG_FILE` environment variable to your config file's path.

```bash
export CONFIG_FILE="/path/to/your/config/file"
```

## Script Overview

The script carries out the following steps:

1. Required libraries are loaded.
2. Project-related variables (topic, name, script path) are defined.
3. Configuration variables are read and the working directory is set.
4. Specific global options are set.
5. The current date is computed.
6. Global functions from external scripts are loaded.
7. The merged table and the HGNC table are loaded for filtering and annotation.

The results are saved and are ready for subsequent analysis or visualization.

### License

This project falls under the MIT license terms. Refer to the LICENSE file for more details.
