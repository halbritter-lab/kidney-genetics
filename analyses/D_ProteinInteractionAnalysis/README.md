
# D. Protein Interaction Analysis

This script (`D_ProteinInteractionAnalysis.R`) is part of the nephrology-focused kidney-genetics project and is designed to analyse interactions between the proteins of the high-evidence kidney-genetic genes (evidence count >= 2).

## Prerequisites

Before running the script, ensure you have the following libraries installed:

- `tidyverse`          (for general table operations)
- `STRINGdb`           (for protein interactions)
- `plotly`             (for generating interactive plots)
- `network`            (for creating protein network objects)
- `ggplot2`            (for creating graphs)
- `R.utils`            (to gzip downloaded and result files)
- `GGally`             (extension to ggplot2)
- `progress`           (for adding progress bars)

## Usage

Run the R script located in the directory where it resides. Adjust paths and configurations based on your specific setup.

```bash
Rscript D_ProteinInteractionAnalysis.R
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
5. Global functions from external scripts are loaded.
6. The current date is computed.
7. A directory for saving results is generated.
8. The annotated HGNC table is loaded.
9. The annotated table of high-evidence genes (evidence count >= 2) is loaded.
10. The table of high-evidence genes is annotated with STRING IDs.
11. The STRING protein interaction database is downloaded and loaded.
12. STRING clusters and subclusters formed by the high-evidence genes are determined and saved.
13. The most probable clinical disease group (tubulopathy, glomerulopathy, cancer, cakut, cyst_cilio, complement, nephrocalcinosis) is extracted for each high-evidence gene, the high-evidence genes are annotated with their respective cluster index and the are results saved.
14. Example plots are generated:
14.1. Plot of clinical disease group distribution within a given subcluster.
14.2. Plot of the interaction network determined by given example index genes, the minimum combined interaction score (STRING) and only direct neighbors.
14.3. Plot of the interaction network determined by given example index genes, the minimum combined interaction score (STRING) and only direct neighbors and first-level indirect neighbors.
14.4. Plot of the interaction network determined by given example index genes, the minimum combined interaction score (STRING) and only direct neighbors.
15. A cluster plot of all high-evidence genes is generated and saved.

### License

This project falls under the MIT license terms. Refer to the LICENSE file for more details.
