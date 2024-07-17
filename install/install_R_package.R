####### Install packages
packExist=require("BiocManager",lib.loc=libLocation,character.only = TRUE)
if (!packExist) {install.packages("BiocManager",repos = "http://cran.rstudio.com/")}
packExist=require("ggplot2",lib.loc=libLocation,character.only = TRUE)
if (!packExist) {install.packages("ggplot2",repos = "http://cran.rstudio.com/")}
BiocManager::install("DNAcopy",force = TRUE)
BiocManager::install("Biostrings",force = TRUE)
########################
