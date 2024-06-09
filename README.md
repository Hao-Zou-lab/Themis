# Themis
![logo](https://github.com/Hao-Zou-lab/Themis/assets/115384996/d133fcff-374a-4710-b832-9712d7231c2f) ![graphical_abstract](https://github.com/Hao-Zou-lab/Themis/assets/115384996/269ad9a2-536f-4fb7-8ad4-3db9788a2755)


Tumor HEterogeneity analysis on Molecular subtypIng System
# Install associated packages
Most dependent packages can be installed by the code below:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

BiocManager::install('NMF')
 
BiocManager::install('GSVA')
 
BiocManager::install('GSEABase')
 
BiocManager::install('ggtree')

install.packages("survival")

install.packages("survminer")

BiocManager::install("limma")

install.packages('dplyr')

install.packages("tidyverse")

install.packages ('sparcl')

install.packages("Rtsne")

install.packages("ggplot2")

install.packages("proxy")

install.packages("pheatmap")

install.packages("ggpubr")

install.packages("reshape2")

install.packages("scales")

install.packages("fmsb")

install.packages("ggprism")

install.packages("devtools")

devtools::install_github ("ricardo-bion/ ggradar ",dependencies = TRUE,force = T")

devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")

# Installation Themis
Installing scSTAR from github: download the file 'Themis_1.0.0.2.tar.gz' and install the package from the local path.

# Preparatory Work
Get the required data ready:
1. The reformatted TCGA and GEO data presented in the paper are available upon request with the proper license.
2. Set the working directory and the "sum data" folder under it. The sum data folder is from the file 'Themis_1.0.0.2.tar.gz' and stores three files: data_meta_sum.RData, tcga_meta_sum.RData, all_pvalue_sum.RData. As shown in the following figure:
   ![figure1](https://github.com/Hao-Zou-lab/Themis/assets/115384996/99f3524a-c9b6-4fc3-9f80-07e77e41bb02)
3. clusterFile: When the user has a custom method, apply it to the directory where the resulting subtyping result file is stored after the platform has existing data. The file naming format is “your method’s name_test data’s name_ cluster.txt”. As shown in the following figure:
   ![figure2](https://github.com/Hao-Zou-lab/Themis/assets/115384996/60dd004b-39d2-4535-b4e4-681b4db224ec)

# Run Themis
library(Themis)

directory <- "~"

path <- "~/sum data"

Themis_custom_method(directory, path)

Themis_process(directory)

Themis_output(directory)

directory is the working directory where the above files (“cluster file”, "sum data" folder, and output results) are stored.
cluster file: When the user has a custom method, apply it to the directory where the resulting subtyping result file is stored after the platform has existing data.

path is the path to the “sum data” folder. The sum data folder stores three files: data_meta_sum.RData, tcga_meta_sum.RData, all_pvalue_sum.RData. data_meta_sum.RData: Summarized the relevant clinical information of the eight literature data mentioned in the article. tcga_meta_sum.RData: Summarized clinical information on 30 different tumor types in TCGA. all_pvalue_sum.RData: The statistical results of the correlation between the subtyping analysis of 38 test data and clinical phenotype using these 8 typing methods. These three files can be obtained from the data folder of the downloaded Themis package.


