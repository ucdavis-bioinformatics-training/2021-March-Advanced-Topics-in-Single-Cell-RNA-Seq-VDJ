### Create a new RStudio project

Open RStudio and create a new project, for more info see (Using-Projects)[https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects]

* File > New Project > New Directory > New Project (name the new directory, Ex. Differential_Expression) and check "use renv with this project" if present.

Learn more about (renv)[https://rstudio.github.io/renv/articles/renv.html]

Set some options and make sure the packages Seurat, sva, ggplot2, dplyr, limma, topGO, WGCNA are installed (if not install it), and then load them and verify they all loaded correctly.

In the R console run the following commands
```r
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if (!any(rownames(installed.packages()) == "rmarkdown")){
  BiocManager::install("rmarkdown")
}

if (!any(rownames(installed.packages()) == "tinytex")){
  BiocManager::install("tinytex")
}

if (!any(rownames(installed.packages()) == "Seurat")){
  BiocManager::install("Seurat")
}

if (!any(rownames(installed.packages()) == "hdf5r")){
  BiocManager::install("hdf5r")
}

if (!any(rownames(installed.packages()) == "knitr")){
  BiocManager::install("knitr")
}

if (!any(rownames(installed.packages()) == "kableExtra")){
  BiocManager::install("kableExtra")
}

if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}

if (!any(rownames(installed.packages()) == "dplyr")){
  BiocManager::install("dplyr")
}

if (!any(rownames(installed.packages()) == "reshape2")){
  BiocManager::install("reshape2")
}

if (!any(rownames(installed.packages()) == "cowplot")){
  BiocManager::install("cowplot")
}

if (!any(rownames(installed.packages()) == "devtools")){
  BiocManager::install("devtools")
}

if (!any(rownames(installed.packages()) == "scRepertoire")){
  devtools::install_github("ncborcherding/scRepertoire")
}

## All of thse should now load without error.

library(rmarkdown)
library(tinytex)
library(Seurat)
library(hdf5r)
library(knitr)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(scRepertoire)

sessionInfo()
```

### Download the template Markdown workshop document VDJ and open it.

In the R console run the following command to download part 1 of data analysis
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-March-Advanced-Topics-in-Single-Cell-RNA-Seq-VDJ/main/data_analysis/VDJ_Analysis.Rmd", "VDJ_Analysis.Rmd")
```

### Download the data for the workshop, extract it.

In the R console run the following command to download and extract the dataset (Little over 1Gb file.

```r
options(timeout=1200)
download.file("https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/4vn7r610cf5d5dv/advsinglecellvdj_March2021.zip", "advsinglecellvdj_March2021.zip")
system("unzip advsinglecellvdj_March2021.zip") # works in Linux and Mac, not sure about Windows"
```

If you timed out on the download, increase 1200 to something higher. If the system command didn't work to extract the zip file, navigate to the folder you downloaded the data in and manually unzip the archive file

**The Dataset will only be available for download during this course**

### Edit the file YAML portion

The top YAML (YAML ain't markup language) portion of the doc tells RStudio how to parse the document.

<pre><code>---
title: "Single Cell VDJ Analysis"
author: your_name
date: current_date
output:
    html_notebook: default
    html_document: default
---</code></pre>


Your RStudio should look something like this

<img src="figures/RStudio.png" alt="RStudio" width="80%"/>


Now spend a few minutes navigating through our data, how may samples are there? Find the hdf5 file and the matrix files. View the html files and lets discuss.
