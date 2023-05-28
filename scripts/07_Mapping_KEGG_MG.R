#####################################################################
## PROJECT : GamfoCyc                                              ##
## STUDIES : APPROVE                                               ##
## AUTHOR : Natacha Koenig                                         ##
## DATE : July 2022                                                ##
## SCRIPT : Récupération des ECs des protéines DE pour les mapper  ##
###         sur les cartes KEGG                                    ##
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
list.of.packages <- c("plyr", "dplyr", "ggplot2", "grid", "gridExtra", "mixOmics", "minfi", "lumi", "stats", "limma", "edgeR", "Heatplus", "made4", "RColorBrewer", "Biobase", "stringr", "readxl", "tibble", "stringr", "tidyr", "pheatmap", "openxlsx") #list of packages required
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] #list of packages non installed
if(length(new.packages)) install.packages(new.packages, repos='https://cran.rstudio.com/') #install packages if the new.packages list is not empty

#tools
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

#PCA and clustering
library(mixOmics)
library(dplyr)
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
library(pheatmap)
library(openxlsx)

#-------------------------------------------------------------------
#  DIRECTORIES                          
#-------------------------------------------------------------------
# The script requires several sub-folders of the current directory :
# /data, /plot and /output

## Working directory
wdir <- getwd()
wdir #current directory
dir()

## Input directories
datadir <- file.path(wdir, "data")

## Output directory
plotdir <- file.path(wdir, "plot")
plotdirACP <- file.path(plotdir, "acp")
plotdirDE <- file.path(plotdir, "DE")

outputdir <- file.path(wdir, "output")
outputdirDE <- file.path(outputdir, "DE")
outputdirKEGG <- file.path(outputdir, "Mapping_KEGG")


##############################################################################
#               I. ECs FROM PROTEOMICS (DE PROTEINS IN ORGANS)               #
##############################################################################
#-------------------------------------------------------------------
#  01. DATA IMPORTATION                           
#-------------------------------------------------------------------
## LOAD DATA FROM DESCRIPTIVE ANALYSIS
load(file=file.path(outputdir, "data.Rdata"), verbose =TRUE)
#all proteins annotated ing gamfocyc and validated by proteomics (n=858)
filt_prot_gamfo_MG_unique 
length(filt_prot_gamfo_MG_unique) # n=858

## LOAD DATA FROM DE ANALYSIS
load(file=file.path(outputdirDE, "data_DE.Rdata"), verbose =TRUE)

#male gonads
resume_gonads_1VS1
DE.down_gonads_1VS1_MG #results DE analysis (script 4)
male_gonads <- DE.down_gonads_1VS1_MG$genes #DE protein list 

#female gonads
resume_gonads_1VS1
DE.up_gonads_1VS1_MG #results DE analysis (script 4)
female_gonads <- DE.up_gonads_1VS1_MG$genes #DE protein list 

#gills
resume_gills_1VSALL
DE.down_gills_1VSALL_MG
gills <- DE.down_gills_1VSALL_MG$genes #DE protein list 

#caeca
resume_caeca_1VSALL
DE.down_caeca_1VSALL_MG
caeca <- DE.down_caeca_1VSALL_MG$genes #DE protein list 


##############################################################################
#               II. ECs FROM TRANSCRIPTOMICS (GAMFOCYC)                      #
##############################################################################
#-------------------------------------------------------------------
#  01. DATA IMPORTATION                           
#-------------------------------------------------------------------
### Import the EC - proteins file from GamfoCyc (fromm all proteins even those that are annotated outside pathways)
raw_EC_proteo <- read.delim(file.path(datadir, "gamfo_functionsbyline.annot"), 
                         check.names=FALSE, stringsAsFactors=FALSE, header=FALSE, skip=1) #lecture table ECs-proteins
raw_EC_proteo <- as_tibble(raw_EC_proteo) #conversion en tibble

raw_EC_proteo <- raw_EC_proteo %>% dplyr::select(V1, V3) #we keep protein name and ECs columns only

# ECs of all proteins annotated in gamfocyc (even outside pathways)
EC_proteo <- raw_EC_proteo %>% rename(V1="Protein_name", V3="EC_number") #rename of the column names

length(unique(EC_proteo$Protein_name)) #62667 : Toutes les protéines même celles qui ne sont pas annotées dans le MG (dans les pathways)


### Filter ECs list from the ECs of proteins that are annotated only in pathways (MG) on gamfocyc 
gamfo_prot_list_unique_MG #liste protéines annotées dans des pathways sur gamfocyc
gamfo_prot_list_unique_MG_clean <- gsub("-PA", "", gamfo_prot_list_unique_MG) #suppression "-PA" dans gamfo_prot_list_unique_MG

EC_prot_in_pathway_gamfocyc <- EC_proteo %>% filter(EC_proteo$Protein_name %in% gamfo_prot_list_unique_MG_clean) #filtre pour récupérer ECs des 4033 protéines dans les pathways gamfocyc



##############################################################################
#       III. RETRIEVE ECs CORRESPONDING TO PROTEINS IN MG (ALL ORGANS)       #
##############################################################################
#-------------------------------------------------------------------
#  01.  ECs FROM TRANSCRIPTOMIC ANNOTATION (GAMFOCYC)                          
#-------------------------------------------------------------------
# col_EC_transcriptomics="#7f7fff"
col_EC_transcriptomics="#A680F3"

list_EC_tr <- str_split(EC_prot_in_pathway_gamfocyc$EC_number, ";")
list_EC_tr <- unlist(list_EC_tr) #69149

trancriptomics_KEGG <- data.frame(EC=list_EC_tr, color=rep(col_EC_transcriptomics, length(list_EC_tr)))
trancriptomics_KEGG <- trancriptomics_KEGG %>% filter(EC!="") #15899

# export table of EC proteins of GM
write.table(trancriptomics_KEGG, file=(file.path(outputdirKEGG, "01_List_EC_all_prot_from_transcriptomics.csv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


#-------------------------------------------------------------------
#  02.  ECs OF PROTEINS ANNOTATED IN PATHWAY (IN MG)                          
#-------------------------------------------------------------------
EC_prot_in_pathway_gamfocyc # 4033 // ECs des protéines annotées à partir données transcriptomiques sur gamfocyc
filt_prot_gamfo_MG_unique # 858 // Protéines annotées à partir de la protéomique et dans le MG 
# col_EC_all_MG="#bfbf7f"
col_EC_all_MG="#19B984"

# intersection of the two tables 
EC_prot_in_MG <- EC_prot_in_pathway_gamfocyc %>% filter(EC_prot_in_pathway_gamfocyc$Protein_name %in% filt_prot_gamfo_MG_unique) #voir pourquoi il manque une protéine 857/858 ? > setdiff(filt_prot_gamfo_MG_unique, EC_prot_in_MG$Protein_name)

list_EC_all_MG <- str_split(EC_prot_in_MG$EC_number, ";")
list_EC_all_MG <- unlist(list_EC_all_MG) #1715

all_prot_MG_KEGG <- data.frame(EC=list_EC_all_MG, color=rep(col_EC_all_MG, length(list_EC_all_MG)))
all_prot_MG_KEGG <- all_prot_MG_KEGG %>% filter(EC!="") #1679

# export table of EC proteins of GM
write.table(all_prot_MG_KEGG, file=(file.path(outputdirKEGG, "01_List_EC_all_prot_MG_from_proteomics.csv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
