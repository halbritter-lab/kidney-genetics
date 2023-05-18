# Kidney-Genetics - Designing a reproducible and curated database of kidney-related genes
Welcome to the GitHub repository  "Kidney-Genetics", a systematically curated, reproducible list of all relevant kidney-related genes known to date. This project aims to provide a unified and standardized database of kidney disease-associated genes, contributing to improved diagnosis, treatment selection, and monitoring of kidney diseases. The database is designed to be automatically updated on a regular basis, ensuring the incorporation of the most up-to-date genetic findings in kidney research.

## Table of contents

- [Overview and Methods](#overview-and-methods)
- [Usage](#usage)
- [File structure](#file-structure)
- [License](#license)
- [Creators](#creators-and-contributors)
- [Contact](#contact)

## Overview and Methods

The Kidney-Genetics Atlas contains information for about 3.000 kidney-associated genes. This information was gathered from various reable sources, such as Genomics England PanelApp, PanelApp Australia, PubTator, OMIM, Orphanet, clinical diagnostic panels, and comprehensive literature review.

## Usage

The Kidney-Genetics Atlas is publicly available and accessible through GitHub. Furthermore, the database is automatically and regularly updated to ensure its currency and relevance.

## File Structure

The repository has the following structure:

```
.
├── analyses/
│   ├── 01_PanelApp/
│   │   ├── data/
│   │   ├── results/
│   │   └── 01_PanelApp.R
│   ├── 02_Literature/
│   │   ├── data/
│   │   ├── results/
│   │   └── 02_Literature.R
│   ├── 03_DiagnosticPanels/
│   │   ├── data/
│   │   ├── results/
│   │   └── 03_DiagnosticPanels.R
│   ├── 04_HPO/
│   │   ├── data/
│   │   ├── results/
│   │   └── 04_HPO.R
│   └── 05_PubTator/
│       ├── data/
│       ├── results/
│       └── 05_PubTator.R
└── functions/
    ├── blueprintgenetics-functions.R
    ├── hgnc-functions.R
    ├── hpo-functions.R
    ├── natera-functions.R
    ├── NCBI-datasets-v2-API-functions.R
    ├── phantomjs-functions.R
    └── PubTator-functions.R
```

- The `analyses/` directory contains the R scripts for different analyses.
- The `functions/` directory contains the necessary functions for HGNC processing.
- The `data/` sub-directory in each analysis folder stores the input data files, including the publication-specific files and the curated overview Excel table.
- The `results/` sub-directory in each analysis folder stores the generated results.

## License

This project is licensed under the terms of the MIT license. For more information, please refer to the [License](LICENSE.md) file.

## Creators and contributors

**Bernt Popp**

- <https://twitter.com/berntpopp>
- <https://github.com/berntpopp>
- <https://orcid.org/0000-0002-3679-1081>
- <https://scholar.google.com/citations?user=Uvhu3t0AAAAJ>

**Constantin Aaron Wolff**

**Nina Rank**

- <https://orcid.org/0000-0002-5984-4836>

**Jan Halbritter**

- <https://orcid.org/0000-0002-1377-9880>

## Contact

If you have any questions, suggestions, or feedback, please feel free to [contact us](contact.md).