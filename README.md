# snSeq_BF
This repository supplies the codes used in the article "Single-nucleus transcriptome atlas of the basal forebrain reveals diverse ageing-related pathways".

If you use this data in your article, please cite: <https://doi.org/10.1093/brain/awaf060>

# Overview


# Scripts Download
The scripts and demo datasets can be downloaded through git:
git clone https://github.com/Scolico/snSeq_BF.git
The download should only take a few seconds.


# Usage Instructions



## System requirements
The majority of these scripts were developed with R/v4.1 on Linux or Windows10 operating systems. 

## R Package Installation
~~~ R
# All dependencies can be dowloaded from CRAN or BiocManager.
#Seurat: 
install.packages('Seurat')
#monocle: 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("monocle")
#GeneSwitches:
list.of.packages <- c("SingleCellExperiment", "Biobase", "fastglm", "ggplot2", "monocle",
                      "plyr", "RColorBrewer", "ggrepel", "ggridges", "gridExtra", "devtools",
                      "mixtools")

## for package "fastglm", "ggplot2", "plyr", "RColorBrewer", "ggrepel", "ggridges", "gridExtra", "mixtools"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## for package "SingleCellExperiment", "Biobase"
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
#Scrublet:
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
#SCENIC:
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
#clusterProfiler /GSVA/hdWGCNA/CellChat/ggplot2
BiocManager::install(c("clusterProfiler","GSVA","hdWGCNA", "CellChat", "ggplot2"))
~~~

# Demo
