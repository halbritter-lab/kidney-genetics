# A. Merge analyses sources

This script (`A_MergeAnalysesSources.R`) is designed to merge and transform analysis sources within the kidney-genetics project. It combines various analysis files and prepares the data for subsequent processing steps.

## Prerequisites

Before running the script, make sure you have the following libraries installed:

- `tidyverse`
- `readr`
- `tools`
- `R.utils`
- `config`

## Usage

Run the R script located in the directory where it resides. You may need to adjust paths and configurations according to your specific setup.

\```bash
Rscript A_MergeAnalysesSources.R
\```

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
5. Merges various analysis files, including defining paths and processing results.

The script saves the merged results in a designated location, ready for further analysis or visualization.

### License

This project is licensed under the terms of the MIT license. See the LICENSE file for details.
