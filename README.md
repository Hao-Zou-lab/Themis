# Themis
![logo](https://github.com/Hao-Zou-lab/Themis/assets/115384996/d133fcff-374a-4710-b832-9712d7231c2f)

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
install.packages("limma") 
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

install.packages("devtools")
devtools::install_github ("ricardo-bion/ ggradar ",dependencies = TRUE,force = T")
devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")

# Installation Themis
install.packages("Themis_1.0.0.1.tar.gz", repos = NULL, type = "source")

# Preparatory Work
Download the required files from data in the current interface:
1. Download the test data file from the pan.baidu.com website (https://pan.baidu.com/s/1dzFwCWckp6Ml6dgkqaJ1tA?pwd=g8le) or the drive.google.com website (https://drive.google.com/file/d/1PITMtvaEsaQFT10vROtQ4uHApXimvCWe/view?usp=sharing).
2. Set the working directory and create a "sum data" folder under it. The sum data folder stores three files: data_meta_sum.RData, tcga_meta_sum.RData, all_pvalue_sum.RData. As shown in the following figure:
   ![figure1](https://github.com/Hao-Zou-lab/Themis/assets/115384996/99f3524a-c9b6-4fc3-9f80-07e77e41bb02)
3. clusterFile: When the user has a custom method, apply it to the directory where the resulting subtyping result file is stored after the platform has existing data. The file naming format is “your method’s name_test data’s name_ cluster.txt”. As shown in the following figure:
   ![figure2](https://github.com/Hao-Zou-lab/Themis/assets/115384996/60dd004b-39d2-4535-b4e4-681b4db224ec)

# Run Themis
Themis_custom_method(directory, path)

Themis_process1(directory)

Themis_output(directory)

directory is the working directory where the above files (“cliFile”, “clusterFile”, “CaculateScore.R”, “Compare_Method.R”, etc.) are stored.
cliFile: Clinical files’ path. The merge results of various analyses (Compare_Method.R) on the sample data, including Stemness Score Caculation (mRNAsi), immune infiltration scores (ESTIMATEScore, ImmuneScore, and StromalScore), immune therapy response prediction (TIDE), as well as clinical information (OS/PFS) and merged clinical risk factors.
clusterFile: When the user has a custom method, apply it to the directory where the resulting subtyping result file is stored after the platform has existing data.

path is the path to the “sum data” folder. The sum data folder stores three files: data_meta_sum.RData, tcga_meta_sum.RData, all_pvalue_sum.RData. data_meta_sum.RData: Summarized the relevant clinical information of the eight literature data mentioned in the article. tcga_meta_sum.RData: Summarized clinical information on 30 different tumor types in TCGA. all_pvalue_sum.RData: The statistical results of the correlation between the subtyping analysis of 38 test data and clinical phenotype using these 8 typing methods. These three files can be obtained from the data folder of the downloaded Themis package.


