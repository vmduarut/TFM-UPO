source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler)
library(scatterplot3d)

source("http://zzlab.net/GAPIT/emma.txt")
setwd("C:/Users/Bisite/Desktop/tfmUPO/GAPIT")

myY <- read.table("C:/Users/Bisite/Desktop/tfmUPO/phenotype1000g.txt", head = TRUE, sep = "\t")
myG <- read.csv("C:/Users/Bisite/Desktop/tfmUPO/ALL.chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.split.AF_0.05.hmp.txt" , head = FALSE, sep = "\t")
myGAPIT_MLMM <- GAPIT(Y=myY, G=myG, model= c("GLM", "MLM", "CMLM", "FarmCPU", "Blink", "FastCPU"), Multiple_analysis=TRUE)
