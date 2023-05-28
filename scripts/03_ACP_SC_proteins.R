#####################################################################
## PROJECT : GamfoCyc                                              ##
## STUDIES : APPROVE                                               ##
## AUTHOR : Amélie Lafont                                          ##
## DATE : July 2022                                                ##
## SCRIPT : Stats descriptives  - Exploration donn?es              ##
#####################################################################

#-------------------------------------------------------------------
#  INDEX.SUMMARY OF THE PIPELINE                   
#-------------------------------------------------------------------

#      PACKAGES & SETTINGS 
#      DIRECTORIES  


#      APPENDICES

#-------------------------------------------------------------------
#  PACKAGES & SETTINGS                          
#-------------------------------------------------------------------
## Required packages
# Installs missing libraries !
list.of.packages <- c("plyr", "dplyr", "ggplot2", "grid", "gridExtra", "mixOmics", "minfi", "lumi", "stats", "limma", "edgeR", "Heatplus", "made4", "RColorBrewer", "Biobase", "stringr", "readxl", "tibble", "stringr", "tidyr") #list of packages required
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] #list of packages non installed
if(length(new.packages)) install.packages(new.packages, repos='https://cran.rstudio.com/') #install packages if the new.packages list is not empty

#tools
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

#PCA and clustering
library(mixOmics)
# library(minfi)
library(dplyr)
# library(lumi)
library(stats)

# Generics
library(limma)
library(edgeR)
library(Heatplus)
library(made4)
library(RColorBrewer)
library(Biobase)
library(stringr)
library(readxl)
library(tibble)
library(stringr)
library(tidyr)
library(lumi)

#-------------------------------------------------------------------
#  DIRECTORIES                          
#-------------------------------------------------------------------
# The script requires several sub-folders of the current directory :
# /data, /plot, /output and /img


## Working directory
wdir <- getwd()
wdir #current directory
dir()

## Input directories
datadir <- file.path(wdir, "data")

## Output directory
plotdir <- file.path(wdir, "plot")
plotdirACP <- file.path(plotdir, "acp")

outputdir <- file.path(wdir, "output")
outputdirACP <- file.path(outputdir, "acp")



###############################################################################
######################### I. ANALYSE POUR ALL PROT GFBM/GFBF ###############
############################## PROTEOCOUNT NON AMBIGUOUS  ####################
###############################################################################
#-------------------------------------------------------------------
#  01. DATA IMPORTATION                           
#-------------------------------------------------------------------
load(file=file.path(outputdir,"data.Rdata"), verbose =TRUE)

# Taking a look to the raw data is always a good idea
names(s2_gfb_counts)
row.names(s2_gfb_counts)
View(s2_gfb_counts)
# read the phenotype data (descriptors, i.e.: exposure, class, batch effects)
colnames(pData)
row.names(pData)<-pData$SampleID
pData<-pData[,-1]

#-------------------------------------------------------------------
#  02. FILTERING LOW ABUNDANCE PROTEINS OR MISSING VALUES      
#                           AND PRE-PROCESSING                      
#-------------------------------------------------------------------
## PRE-FILTERING STEP
# Put the data into a DGEList object (edgeR) y
y <- DGEList(counts=s2_gfb_counts, genes=rownames(s2_gfb_counts))
dim(y)   #[1] 4150   36
n_prot <- dim(y)[[1]]
y$samples
y$genes

## LOOK AT THE DATA
## Euclidian clustering
tiff(file=file.path(plotdirACP, "01_CAH_before_filtering.tiff"))
plotSampleRelation(y$counts, method = "cluster", labels=pData$SexOrgan, cex=1, main=paste("Hierarchical clustering based on", n_prot,"proteins raw data"))
dev.off()

## PCA
ty<-as.data.frame(t(y$counts)) #transpose count dataframe y

# Exploring how many dimensions explain variability
tiff(file=file.path(plotdirACP, "02_PCA_axes_before_filtering.tiff"))
plot(tune.pca(ty, ncomp=10, center=TRUE, scale=FALSE), main="PCA axes before filtering") #explore how many dimensions explain variability #here 2
dev.off()

pca =pca(ty, ncomp=2)

# plotting the PCA results
tiff(file=file.path(plotdirACP, "03_PCA_before_filtering.tiff"))
plotIndiv(pca,
          group = pData$SexOrgan,
          col.per.group = c("#32B866", "#8FE0AE", "#3498DB", "#96CAED", "#A624DD", "#C573E8", "#CB001B", "#F6798A", "#FF9317", "#FFBC70", "#3A2865", "#543A92"),
          ind.names = F,
          legend=T, legend.position="right", legend.title=NULL,
          size.xlabel = 15,
          size.ylabel = 15,
          size.legend = 15,
          size.axis = 15,
          size.title = 15,
          title='Samples PCA (Organ)',
          comp=c(1,2),
          ellipse=F,
          pch=16)

dev.off()

## FILTERING
# Define a threshold corresponding to a count of 3 spectral counts for filtering
sc <- 3
th <- cpm(sc, mean(colSums(s2_gfb_counts)))[1] #converting SC to cpm
round(th) #[1] 134
nbsamples <- 3 #mimimum number of samples threshold
keep <- rowSums(cpm(s2_gfb_counts) > th) >=nbsamples #filtering step
s2_gfb_countsF <- s2_gfb_counts[keep,]
dim(s2_gfb_countsF)
#[1] 2190   36

#-------------------------------------------------------------------
#  03. DATA NORMALIZATION                      
#-------------------------------------------------------------------
# Normalization by the calcNormfactors function of the edgeR package
# based on a TMM (Trimmed mean of M-values) normalization procedure
#(https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)

group = as.factor(pData$SexOrgan)

# Put the filtered data into a DGEList object yF
yF<-DGEList(counts=s2_gfb_countsF,
            genes=row.names(s2_gfb_countsF),
            group=group)
dim(yF)
head(yF$samples) #containing sample informations

# Normalization (required for WGCNA) yN
yN <- calcNormFactors(yF) #normalization of the filtered counts 
dim(yN)
head(yN$samples) #containing sample informations

# Create the data with counts multiplied by yN norm.factors : yN2
yN2 <- yN
n <- length(colnames(yN$counts)) #number of columns of the yN counts table
n #number of samples

for (i in seq (from=1, to=n)) {
  yN2$counts[,i] <- yN$counts[,i]*yN$samples$norm.factors[i]
}

dim(yN2)
nF_prot <- dim(yN2)[[1]]
head(yN2$genes)
data <- yN2$counts #normalized counts data yN2 for the rest of the analysis
rownames(data)
head(data)

#-------------------------------------------------------------------
#  04.  EXPLORATORY DATA ANALYSIS                   
#-------------------------------------------------------------------
##  LOOK AT THE DATA AFTER FILTERING AND NORMALIZATION
## Euclidian clustering
tiff(file=file.path(plotdirACP, "04_CAH_after_filt_norm.tiff"))
plotSampleRelation(yN2$counts, method = "cluster",labels=pData$SexOrgan, cex=1, main=paste("Hierarchical clustering based on", nF_prot, "proteins \n filtered and normalized data")) #no more outliers
dev.off()

## PCA
tyN2<-as.data.frame(t(yN2$counts)) #transpose count dataframe yN2

# Exploring how many dimensions explain variability
tiff(file=file.path(plotdirACP, "05_PCA_axes_after_filt_norm.tiff"))
plot(tune.pca(tyN2, ncomp=10, center=TRUE, scale=FALSE), main="PCA axes after filtering") #here 2 dimensions
dev.off()

pca =pca(tyN2, ncomp=2) #we choose 2 dimensions for the PCA

# plotting the PCA results
tiff(file=file.path(plotdirACP, "06_PCA_after_filt_norm.tiff"))
plotIndiv(pca,
          group = pData$SexOrgan,
          col.per.group = c("#32B866", "#8FE0AE", "#3498DB", "#96CAED", "#A624DD", "#C573E8", "#CB001B", "#F6798A", "#FF9317", "#FFBC70", "#3A2865", "#543A92"),
          ind.names = F,
          legend=T, legend.position="right", legend.title=NULL,
          size.xlabel = 15,
          size.ylabel = 15,
          size.legend = 15,
          size.axis = 15,
          size.title = 15,
          title='PCA filtered and normalized samples',
          comp=c(1,2),
          ellipse=F,
          pch=16)
dev.off()


#-------------------------------------------------------------------
#  05.  SAVING DATA FOR SCRIPT N.04
#-------------------------------------------------------------------
save(s2_gfb_countsF, pData, yN, group, yN2, file=file.path(outputdirACP,"data_ACP.Rdata"))

