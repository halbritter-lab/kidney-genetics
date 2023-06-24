# Kidney-Genetics Documentation

The repository subfolder for the Kidney-Genetics documentation.

Tto build the documentation execute following commands:

```
## load libraries
library(bookdown)
library(config)

project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/edit_docs/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

# render the book
bookdown::render_book("index.Rmd", "all")
```

Rendering HTML widgets for PDF requires webshot and phantomjs [FROM:](https://bookdown.org/yihui/bookdown/html-widgets.html).
```
install.packages("webshot")
webshot::install_phantomjs()
```

The CSL file is [American Psychological Association 7th edition, from Zotero](https://www.zotero.org/styles/apa).

## TODO
- automatic loading and filtering of current result tables in all Rmd files
- change scripts to load the gzipped files
- make a script that runs the above commands to re-generate the pages website