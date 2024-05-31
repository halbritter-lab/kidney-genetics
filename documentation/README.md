# Kidney-Genetics Documentation Markdown and Scripts

The repository subfolder for the Kidney-Genetics documentation.

<code>The documentation is automatically built and deployed to https://halbritter-lab.github.io/kidney-genetics/ using GitHub Actions into the branch hg-pages.
The configuration for the GitHub Pages build is in the `.github/workflows` folder.</code>

For local development only, the documentation can be built using the following commands:

```
## load libraries
library(bookdown)
library(config)

project_topic <- "nephrology"
project_name <- "kidney-genetics"
script_path <- "/documentation/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

# render the book
bookdown::render_book("index.Rmd", "all")
```


The CSL file is [American Psychological Association 7th edition, from Zotero](https://www.zotero.org/styles/apa).

## required packages
Implementing the draw.io files requires the [knitrdrawio package](https://github.com/rchaput/knitrdrawio).
See installation instructions there or use following commands:

```
# Install the remotes package to install packages directly from GitHub
install.packages("remotes")
# Now, install the knitrdrawio package from GitHub
remotes::install_github("rchaput/knitrdrawio")
```

The knitrdrawio package also requires installing [draw.io](https://github.com/jgraph/drawio-desktop/releases) and setting the path to the executable in the YAML config file.

Rendering HTML widgets for PDF requires webshot and phantomjs [FROM:](https://bookdown.org/yihui/bookdown/html-widgets.html).
```
install.packages("webshot")
webshot::install_phantomjs()
```