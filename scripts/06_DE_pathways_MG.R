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
outputdirACP <- file.path(outputdir, "acp")
outputdirDE <- file.path(outputdir, "DE")
outputdirPATH <- file.path(outputdir, "Pathways/MG")
outputdirPATH_gills <- file.path(outputdirPATH, "gills")
outputdirPATH_caeca <- file.path(outputdirPATH, "caeca")
outputdirPATH_gonads <- file.path(outputdirPATH, "gonads")

#-------------------------------------------------------------------
#  LOAD DATA                          
#-------------------------------------------------------------------
load(file=file.path(outputdir, "data.Rdata"), verbose =TRUE)
load(file=file.path(outputdirDE, "data_DE.Rdata"), verbose =TRUE)
resume_gills_1VS1
resume_gills_1VSALL
resume_caeca_1VS1
resume_caeca_1VSALL
resume_gonads_1VS1
resume_gonads_1VSALL

###############################################################################
################################ PATHWAY PROT DE###############################
###############################################################################

# #############################################################################
#                 I. ANALYSE DIFFERENTIELLE GFB 1 vs 1 ORGANS                 #
# #############################################################################

#################################### GILLS ####################################  

#-------------------------------------------------------------------
#  LOAD DATA : LIST MG PROT DE                          
#-------------------------------------------------------------------
# ALL PROT
all_prot_DE_gills_1VS1$genes <- str_extract_all(all_prot_DE_gills_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gills_1VS1 <- unlist(all_prot_DE_gills_1VS1$genes)
length(protDE_gills_1VS1)

## UP
DE.up_gills_1VS1$genes <- str_extract_all(DE.up_gills_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gills_1VS1 <- unlist(DE.up_gills_1VS1$genes)
length(UP_list_gills_1VS1)

## DOWN
DE.down_gills_1VS1$genes <- str_extract_all(DE.down_gills_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gills_1VS1 <- unlist(DE.down_gills_1VS1$genes)
length(DOWN_list_gills_1VS1)


################ MG ################  
# ALL PROT
all_prot_DE_gills_1VS1_MG$genes <- str_extract_all(all_prot_DE_gills_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gills_1VS1_MG <- unlist(all_prot_DE_gills_1VS1_MG$genes)
length(protDE_gills_1VS1_MG)

## UP
DE.up_gills_1VS1_MG$genes <- str_extract_all(DE.up_gills_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gills_1VS1_MG <- unlist(DE.up_gills_1VS1_MG$genes)
length(UP_list_gills_1VS1_MG)

## DOWN
DE.down_gills_1VS1_MG$genes <- str_extract_all(DE.down_gills_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gills_1VS1_MG <- unlist(DE.down_gills_1VS1_MG$genes)
length(DOWN_list_gills_1VS1_MG)



#----------------------------------------------------------------------------------
#  RETRIEVE PATHWAY OF MG PROT INVOLVED IN + OCCURENCE/FREQUENCE PATHWAY + EXPORT                       
#----------------------------------------------------------------------------------
#ALL DE PROT
prot_DE_pathway_gills_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gills_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gills_1VS1_all <- prot_DE_pathway_gills_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gills_1VS1_all, file = file.path(outputdirPATH_gills, "04_prot_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gills_1VS1_all, file = file.path(outputdirPATH_gills,"04b_group_ALL_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gills_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gills_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gills_1VS1_all <- prot_DOWN_pathway_gills_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gills_1VS1_all, file = file.path(outputdirPATH_gills, "05_prot_DOWN_pathway_gills_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gills_1VS1_all, file = file.path(outputdirPATH_gills, "05b_group_DOWN_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gills_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gills_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gills_1VS1_all <- prot_UP_pathway_gills_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gills_1VS1_all, file = file.path(outputdirPATH_gills, "06_prot_UP_pathway_gills_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gills_1VS1_all, file = file.path(outputdirPATH_gills, "06b_group_UP_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)


################ MG ################  
#ALL DE PROT
prot_DE_pathway_gills_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gills_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gills_1VS1_MG <- prot_DE_pathway_gills_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gills_1VS1_MG, file = file.path(outputdirPATH_gills, "04_prot_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gills_1VS1_MG, file = file.path(outputdirPATH_gills,"04b_group_ALL_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gills_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gills_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gills_1VS1_MG <- prot_DOWN_pathway_gills_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gills_1VS1_MG, file = file.path(outputdirPATH_gills, "05_prot_DOWN_pathway_gills_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gills_1VS1_MG, file = file.path(outputdirPATH_gills, "05b_group_DOWN_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gills_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gills_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gills_1VS1_MG <- prot_UP_pathway_gills_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gills_1VS1_MG, file = file.path(outputdirPATH_gills, "06_prot_UP_pathway_gills_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gills_1VS1_MG, file = file.path(outputdirPATH_gills, "06b_group_UP_DE_pathway_gills_1VS1.xlsx"), overwrite=TRUE)


#################################### CAECA ####################################  

#-------------------------------------------------------------------
#  LOAD DATA : LIST MG PROT DE                          
#-------------------------------------------------------------------
# ALL PROT
all_prot_DE_caeca_1VS1$genes <- str_extract_all(all_prot_DE_caeca_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_caeca_1VS1 <- unlist(all_prot_DE_caeca_1VS1$genes)
length(protDE_caeca_1VS1)

## UP
DE.up_caeca_1VS1$genes <- str_extract_all(DE.up_caeca_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_caeca_1VS1 <- unlist(DE.up_caeca_1VS1$genes)
length(UP_list_caeca_1VS1)

## DOWN
DE.down_caeca_1VS1$genes <- str_extract_all(DE.down_caeca_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_caeca_1VS1 <- unlist(DE.down_caeca_1VS1$genes)
length(DOWN_list_caeca_1VS1)

################ MG ################  
# ALL PROT
all_prot_DE_caeca_1VS1_MG$genes <- str_extract_all(all_prot_DE_caeca_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_caeca_1VS1_MG <- unlist(all_prot_DE_caeca_1VS1_MG$genes)
length(protDE_caeca_1VS1_MG)

## UP
DE.up_caeca_1VS1_MG$genes <- str_extract_all(DE.up_caeca_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_caeca_1VS1_MG <- unlist(DE.up_caeca_1VS1_MG$genes)
length(UP_list_caeca_1VS1_MG)

## DOWN
DE.down_caeca_1VS1_MG$genes <- str_extract_all(DE.down_caeca_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_caeca_1VS1_MG <- unlist(DE.down_caeca_1VS1_MG$genes)
length(DOWN_list_caeca_1VS1_MG)

#----------------------------------------------------------------------------------
#  RETRIEVE PATHWAY OF MG PROT INVOLVED IN + OCCURENCE/FREQUENCE PATHWAY + EXPORT                       
#----------------------------------------------------------------------------------
#ALL DE PROT
prot_DE_pathway_caeca_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_caeca_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_caeca_1VS1_all <- prot_DE_pathway_caeca_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_caeca_1VS1_all, file = file.path(outputdirPATH_caeca, "04_prot_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_caeca_1VS1_all, file = file.path(outputdirPATH_caeca,"04b_group_ALL_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_caeca_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_caeca_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_caeca_1VS1_all <- prot_DOWN_pathway_caeca_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_caeca_1VS1_all, file = file.path(outputdirPATH_caeca, "05_prot_DOWN_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_caeca_1VS1_all, file = file.path(outputdirPATH_caeca, "05b_group_DOWN_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_caeca_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_caeca_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_caeca_1VS1_all <- prot_UP_pathway_caeca_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_caeca_1VS1_all, file = file.path(outputdirPATH_caeca, "06_prot_UP_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_caeca_1VS1_all, file = file.path(outputdirPATH_caeca, "06b_group_UP_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)

################ MG ################  
#ALL DE PROT
prot_DE_pathway_caeca_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_caeca_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_caeca_1VS1_MG <- prot_DE_pathway_caeca_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_caeca_1VS1_MG, file = file.path(outputdirPATH_caeca, "04_prot_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_caeca_1VS1_MG, file = file.path(outputdirPATH_caeca,"04b_group_ALL_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_caeca_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_caeca_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_caeca_1VS1_MG <- prot_DOWN_pathway_caeca_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_caeca_1VS1_MG, file = file.path(outputdirPATH_caeca, "05_prot_DOWN_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_caeca_1VS1_MG, file = file.path(outputdirPATH_caeca, "05b_group_DOWN_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)


#UP PROT
prot_UP_pathway_caeca_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_caeca_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_caeca_1VS1_MG <- prot_UP_pathway_caeca_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_caeca_1VS1_MG, file = file.path(outputdirPATH_caeca, "06_prot_UP_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_caeca_1VS1_MG, file = file.path(outputdirPATH_caeca, "06b_group_UP_DE_pathway_caeca_1VS1.xlsx"), overwrite=TRUE)

#################################### GONADES ####################################  

#-------------------------------------------------------------------
#  LOAD DATA : LIST MG PROT DE                          
#-------------------------------------------------------------------
# ALL PROT
all_prot_DE_gonads_1VS1$genes <- str_extract_all(all_prot_DE_gonads_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gonads_1VS1 <- unlist(all_prot_DE_gonads_1VS1$genes)
length(protDE_gonads_1VS1)

## UP
DE.up_gonads_1VS1$genes <- str_extract_all(DE.up_gonads_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gonads_1VS1 <- unlist(DE.up_gonads_1VS1$genes)
length(UP_list_gonads_1VS1)

## DOWN
DE.down_gonads_1VS1$genes <- str_extract_all(DE.down_gonads_1VS1$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gonads_1VS1 <- unlist(DE.down_gonads_1VS1$genes)
length(DOWN_list_gonads_1VS1)


################ MG ################  
# ALL PROT
all_prot_DE_gonads_1VS1_MG$genes <- str_extract_all(all_prot_DE_gonads_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gonads_1VS1_MG <- unlist(all_prot_DE_gonads_1VS1_MG$genes)
length(protDE_gonads_1VS1_MG)

## UP
DE.up_gonads_1VS1_MG$genes <- str_extract_all(DE.up_gonads_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gonads_1VS1_MG <- unlist(DE.up_gonads_1VS1_MG$genes)
length(UP_list_gonads_1VS1_MG)

## DOWN
DE.down_gonads_1VS1_MG$genes <- str_extract_all(DE.down_gonads_1VS1_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gonads_1VS1_MG <- unlist(DE.down_gonads_1VS1_MG$genes)
length(DOWN_list_gonads_1VS1_MG)



#----------------------------------------------------------------------------------
#  RETRIEVE PATHWAY OF MG PROT INVOLVED IN + OCCURENCE/FREQUENCE PATHWAY + EXPORT                       
#----------------------------------------------------------------------------------
#ALL DE PROT
prot_DE_pathway_gonads_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gonads_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gonads_1VS1_all <- prot_DE_pathway_gonads_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gonads_1VS1_all, file = file.path(outputdirPATH_gonads, "04_prot_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gonads_1VS1_all, file = file.path(outputdirPATH_gonads,"04b_group_ALL_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gonads_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gonads_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gonads_1VS1_all <- prot_DOWN_pathway_gonads_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gonads_1VS1_all, file = file.path(outputdirPATH_gonads, "05_prot_DOWN_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gonads_1VS1_all, file = file.path(outputdirPATH_gonads, "05b_group_DOWN_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gonads_1VS1_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gonads_1VS1)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gonads_1VS1_all <- prot_UP_pathway_gonads_1VS1_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gonads_1VS1_all, file = file.path(outputdirPATH_gonads, "06_prot_UP_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gonads_1VS1_all, file = file.path(outputdirPATH_gonads, "06b_group_UP_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)


################ MG ################  
#ALL DE PROT
prot_DE_pathway_gonads_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gonads_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gonads_1VS1_MG <- prot_DE_pathway_gonads_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gonads_1VS1_MG, file = file.path(outputdirPATH_gonads, "04_prot_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gonads_1VS1_MG, file = file.path(outputdirPATH_gonads,"04b_group_ALL_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gonads_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gonads_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gonads_1VS1_MG <- prot_DOWN_pathway_gonads_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gonads_1VS1_MG, file = file.path(outputdirPATH_gonads, "05_prot_DOWN_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gonads_1VS1_MG, file = file.path(outputdirPATH_gonads, "05b_group_DOWN_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gonads_1VS1_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gonads_1VS1_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gonads_1VS1_MG <- prot_UP_pathway_gonads_1VS1_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gonads_1VS1_MG, file = file.path(outputdirPATH_gonads, "06_prot_UP_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gonads_1VS1_MG, file = file.path(outputdirPATH_gonads, "06b_group_UP_DE_pathway_gonads_1VS1.xlsx"), overwrite=TRUE)

###############################################################################
#                II. ANALYSE DIFFERENTIELLE GFB 1 vs ALL ORGANS               #
###############################################################################

#################################### GILLS ####################################  

#-------------------------------------------------------------------
#  LOAD DATA : LIST MG PROT DE                          
#-------------------------------------------------------------------
# ALL PROT
all_prot_DE_gills_1VSALL$genes <- str_extract_all(all_prot_DE_gills_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gills_1VSALL <- unlist(all_prot_DE_gills_1VSALL$genes)
length(protDE_gills_1VSALL)

## UP
DE.up_gills_1VSALL$genes <- str_extract_all(DE.up_gills_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gills_1VSALL <- unlist(DE.up_gills_1VSALL$genes)
length(UP_list_gills_1VSALL)

## UP.10
DE.up.10_gills_1VSALL$genes <- str_extract_all(DE.up.10_gills_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE up
UP10_list_gills_1VSALL <- unlist(DE.up.10_gills_1VSALL$genes)
length(UP10_list_gills_1VSALL)

## DOWN
DE.down_gills_1VSALL$genes <- str_extract_all(DE.down_gills_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gills_1VSALL <- unlist(DE.down_gills_1VSALL$genes)
length(DOWN_list_gills_1VSALL)

## DOWN.10
DE.down.10_gills_1VSALL$genes <- str_extract_all(DE.down.10_gills_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE down
DOWN10_list_gills_1VSALL <- unlist(DE.down.10_gills_1VSALL$genes)
length(DOWN10_list_gills_1VSALL)

################ MG ################  
# ALL PROT
all_prot_DE_gills_1VSALL_MG$genes <- str_extract_all(all_prot_DE_gills_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gills_1VSALL_MG <- unlist(all_prot_DE_gills_1VSALL_MG$genes)
length(protDE_gills_1VSALL_MG)

## UP
DE.up_gills_1VSALL_MG$genes <- str_extract_all(DE.up_gills_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gills_1VSALL_MG <- unlist(DE.up_gills_1VSALL_MG$genes)
length(UP_list_gills_1VSALL_MG)

## UP.10
DE.up.10_gills_1VSALL_MG$genes <- str_extract_all(DE.up.10_gills_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE up
UP10_list_gills_1VSALL_MG <- unlist(DE.up.10_gills_1VSALL_MG$genes)
length(UP10_list_gills_1VSALL_MG)

## DOWN
DE.down_gills_1VSALL_MG$genes <- str_extract_all(DE.down_gills_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gills_1VSALL_MG <- unlist(DE.down_gills_1VSALL_MG$genes)
length(DOWN_list_gills_1VSALL_MG)

## DOWN.10
DE.down.10_gills_1VSALL_MG$genes <- str_extract_all(DE.down.10_gills_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE down
DOWN10_list_gills_1VSALL_MG <- unlist(DE.down.10_gills_1VSALL_MG$genes)
length(DOWN10_list_gills_1VSALL_MG)


#----------------------------------------------------------------------------------
#  RETRIEVE PATHWAY OF MG PROT INVOLVED IN + OCCURENCE/FREQUENCE PATHWAY + EXPORT                       
#----------------------------------------------------------------------------------
#ALL DE PROT
prot_DE_pathway_gills_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gills_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gills_1VSALL_all <- prot_DE_pathway_gills_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "04_prot_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gills_1VSALL_all, file = file.path(outputdirPATH_gills,"04b_group_ALL_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gills_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gills_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gills_1VSALL_all <- prot_DOWN_pathway_gills_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "05_prot_DOWN_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "05b_group_DOWN_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

# #DOWN10 PROT
prot_DOWN10_pathway_gills_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN10_list_gills_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN10_gills_1VSALL_all <- prot_DOWN10_pathway_gills_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN10_pathway_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "05_prot_DOWN10_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN10_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "05b_group_DOWN10_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gills_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gills_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gills_1VSALL_all <- prot_UP_pathway_gills_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "06_prot_UP_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "06b_group_UP_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

#UP10 PROT
prot_UP10_pathway_gills_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP10_list_gills_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP10_gills_1VSALL_all <- prot_UP10_pathway_gills_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP10_pathway_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "06_prot_UP10_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP10_gills_1VSALL_all, file = file.path(outputdirPATH_gills, "06b_group_UP10_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

################ MG ################  
#ALL DE PROT
prot_DE_pathway_gills_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gills_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gills_1VSALL_MG <- prot_DE_pathway_gills_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "04_prot_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gills_1VSALL_MG, file = file.path(outputdirPATH_gills,"04b_group_ALL_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gills_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gills_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gills_1VSALL_MG <- prot_DOWN_pathway_gills_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "05_prot_DOWN_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "05b_group_DOWN_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

# #DOWN10 PROT
prot_DOWN10_pathway_gills_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN10_list_gills_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN10_gills_1VSALL_MG <- prot_DOWN10_pathway_gills_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN10_pathway_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "05_prot_DOWN10_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN10_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "05b_group_DOWN10_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gills_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gills_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gills_1VSALL_MG <- prot_UP_pathway_gills_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "06_prot_UP_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "06b_group_UP_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)

#UP10 PROT
prot_UP10_pathway_gills_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP10_list_gills_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP10_gills_1VSALL_MG <- prot_UP10_pathway_gills_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP10_pathway_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "06_prot_UP10_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP10_gills_1VSALL_MG, file = file.path(outputdirPATH_gills, "06b_group_UP10_DE_pathway_gills_1VSALL.xlsx"), overwrite=TRUE)


#################################### CAECA ####################################  

#-------------------------------------------------------------------
#  LOAD DATA : LIST MG PROT DE                          
#-------------------------------------------------------------------
# ALL PROT
all_prot_DE_caeca_1VSALL$genes <- str_extract_all(all_prot_DE_caeca_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_caeca_1VSALL <- unlist(all_prot_DE_caeca_1VSALL$genes)
length(protDE_caeca_1VSALL)

## UP
DE.up_caeca_1VSALL$genes <- str_extract_all(DE.up_caeca_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_caeca_1VSALL <- unlist(DE.up_caeca_1VSALL$genes)
length(UP_list_caeca_1VSALL)

## UP.10
DE.up.10_caeca_1VSALL$genes <- str_extract_all(DE.up.10_caeca_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE up
UP10_list_caeca_1VSALL <- unlist(DE.up.10_caeca_1VSALL$genes)
length(UP10_list_caeca_1VSALL)

## DOWN
DE.down_caeca_1VSALL$genes <- str_extract_all(DE.down_caeca_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_caeca_1VSALL <- unlist(DE.down_caeca_1VSALL$genes)
length(DOWN_list_caeca_1VSALL)

## DOWN.10
DE.down.10_caeca_1VSALL$genes <- str_extract_all(DE.down.10_caeca_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE down
DOWN10_list_caeca_1VSALL <- unlist(DE.down.10_caeca_1VSALL$genes)
length(DOWN10_list_caeca_1VSALL)

################ MG ################  
# ALL PROT
all_prot_DE_caeca_1VSALL_MG$genes <- str_extract_all(all_prot_DE_caeca_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_caeca_1VSALL_MG <- unlist(all_prot_DE_caeca_1VSALL_MG$genes)
length(protDE_caeca_1VSALL_MG)

## UP
DE.up_caeca_1VSALL_MG$genes <- str_extract_all(DE.up_caeca_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_caeca_1VSALL_MG <- unlist(DE.up_caeca_1VSALL_MG$genes)
length(UP_list_caeca_1VSALL_MG)

## UP.10
DE.up.10_caeca_1VSALL_MG$genes <- str_extract_all(DE.up.10_caeca_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE up
UP10_list_caeca_1VSALL_MG <- unlist(DE.up.10_caeca_1VSALL_MG$genes)
length(UP10_list_caeca_1VSALL_MG)

## DOWN
DE.down_caeca_1VSALL_MG$genes <- str_extract_all(DE.down_caeca_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_caeca_1VSALL_MG <- unlist(DE.down_caeca_1VSALL_MG$genes)
length(DOWN_list_caeca_1VSALL_MG)

## DOWN.10
DE.down.10_caeca_1VSALL_MG$genes <- str_extract_all(DE.down.10_caeca_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE down
DOWN10_list_caeca_1VSALL_MG <- unlist(DE.down.10_caeca_1VSALL_MG$genes)
length(DOWN10_list_caeca_1VSALL_MG)


#----------------------------------------------------------------------------------
#  RETRIEVE PATHWAY OF MG PROT INVOLVED IN + OCCURENCE/FREQUENCE PATHWAY + EXPORT                       
#----------------------------------------------------------------------------------
#ALL DE PROT
prot_DE_pathway_caeca_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_caeca_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_caeca_1VSALL_all <- prot_DE_pathway_caeca_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "04_prot_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca,"04b_group_ALL_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_caeca_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_caeca_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_caeca_1VSALL_all <- prot_DOWN_pathway_caeca_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "05_prot_DOWN_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "05b_group_DOWN_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

# #DOWN10 PROT
prot_DOWN10_pathway_caeca_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN10_list_caeca_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN10_caeca_1VSALL_all <- prot_DOWN10_pathway_caeca_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN10_pathway_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "05_prot_DOWN10_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN10_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "05b_group_DOWN10_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_caeca_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_caeca_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_caeca_1VSALL_all <- prot_UP_pathway_caeca_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "06_prot_UP_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "06b_group_UP_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

#UP10 PROT
prot_UP10_pathway_caeca_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP10_list_caeca_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP10_caeca_1VSALL_all <- prot_UP10_pathway_caeca_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP10_pathway_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "06_prot_UP10_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP10_caeca_1VSALL_all, file = file.path(outputdirPATH_caeca, "06b_group_UP10_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

################ MG ################  
#ALL DE PROT
prot_DE_pathway_caeca_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_caeca_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_caeca_1VSALL_MG <- prot_DE_pathway_caeca_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "04_prot_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca,"04b_group_ALL_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_caeca_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_caeca_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_caeca_1VSALL_MG <- prot_DOWN_pathway_caeca_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "05_prot_DOWN_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "05b_group_DOWN_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

# #DOWN10 PROT
prot_DOWN10_pathway_caeca_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN10_list_caeca_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN10_caeca_1VSALL_MG <- prot_DOWN10_pathway_caeca_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN10_pathway_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "05_prot_DOWN10_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN10_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "05b_group_DOWN10_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_caeca_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_caeca_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_caeca_1VSALL_MG <- prot_UP_pathway_caeca_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "06_prot_UP_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "06b_group_UP_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

#UP10 PROT
prot_UP10_pathway_caeca_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP10_list_caeca_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP10_caeca_1VSALL_MG <- prot_UP10_pathway_caeca_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP10_pathway_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "06_prot_UP10_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP10_caeca_1VSALL_MG, file = file.path(outputdirPATH_caeca, "06b_group_UP10_DE_pathway_caeca_1VSALL.xlsx"), overwrite=TRUE)

#################################### GONADS ####################################  

#-------------------------------------------------------------------
#  LOAD DATA : LIST MG PROT DE                          
#-------------------------------------------------------------------
# ALL PROT
all_prot_DE_gonads_1VSALL$genes <- str_extract_all(all_prot_DE_gonads_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gonads_1VSALL <- unlist(all_prot_DE_gonads_1VSALL$genes)
length(protDE_gonads_1VSALL)

## UP
DE.up_gonads_1VSALL$genes <- str_extract_all(DE.up_gonads_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gonads_1VSALL <- unlist(DE.up_gonads_1VSALL$genes)
length(UP_list_gonads_1VSALL)

## UP.10
DE.up.10_gonads_1VSALL$genes <- str_extract_all(DE.up.10_gonads_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE up
UP10_list_gonads_1VSALL <- unlist(DE.up.10_gonads_1VSALL$genes)
length(UP10_list_gonads_1VSALL)

## DOWN
DE.down_gonads_1VSALL$genes <- str_extract_all(DE.down_gonads_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gonads_1VSALL <- unlist(DE.down_gonads_1VSALL$genes)
length(DOWN_list_gonads_1VSALL)

## DOWN.10
DE.down.10_gonads_1VSALL$genes <- str_extract_all(DE.down.10_gonads_1VSALL$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE down
DOWN10_list_gonads_1VSALL <- unlist(DE.down.10_gonads_1VSALL$genes)
length(DOWN10_list_gonads_1VSALL)

################ MG ################  
# ALL PROT
all_prot_DE_gonads_1VSALL_MG$genes <- str_extract_all(all_prot_DE_gonads_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
protDE_gonads_1VSALL_MG <- unlist(all_prot_DE_gonads_1VSALL_MG$genes)
length(protDE_gonads_1VSALL_MG)

## UP
DE.up_gonads_1VSALL_MG$genes <- str_extract_all(DE.up_gonads_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
UP_list_gonads_1VSALL_MG <- unlist(DE.up_gonads_1VSALL_MG$genes)
length(UP_list_gonads_1VSALL_MG)

## UP.10
DE.up.10_gonads_1VSALL_MG$genes <- str_extract_all(DE.up.10_gonads_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE up
UP10_list_gonads_1VSALL_MG <- unlist(DE.up.10_gonads_1VSALL_MG$genes)
length(UP10_list_gonads_1VSALL_MG)

## DOWN
DE.down_gonads_1VSALL_MG$genes <- str_extract_all(DE.down_gonads_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}")
DOWN_list_gonads_1VSALL_MG <- unlist(DE.down_gonads_1VSALL_MG$genes)
length(DOWN_list_gonads_1VSALL_MG)

## DOWN.10
DE.down.10_gonads_1VSALL_MG$genes <- str_extract_all(DE.down.10_gonads_1VSALL_MG$genes, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}") #Liste top 10 proteines DE down
DOWN10_list_gonads_1VSALL_MG <- unlist(DE.down.10_gonads_1VSALL_MG$genes)
length(DOWN10_list_gonads_1VSALL_MG)


#----------------------------------------------------------------------------------
#  RETRIEVE PATHWAY OF MG PROT INVOLVED IN + OCCURENCE/FREQUENCE PATHWAY + EXPORT                       
#----------------------------------------------------------------------------------
#ALL DE PROT
prot_DE_pathway_gonads_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gonads_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gonads_1VSALL_all <- prot_DE_pathway_gonads_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "04_prot_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads,"04b_group_ALL_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gonads_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gonads_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gonads_1VSALL_all <- prot_DOWN_pathway_gonads_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "05_prot_DOWN_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "05b_group_DOWN_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

# #DOWN10 PROT
prot_DOWN10_pathway_gonads_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN10_list_gonads_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN10_gonads_1VSALL_all <- prot_DOWN10_pathway_gonads_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN10_pathway_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "05_prot_DOWN10_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN10_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "05b_group_DOWN10_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gonads_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gonads_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gonads_1VSALL_all <- prot_UP_pathway_gonads_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "06_prot_UP_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "06b_group_UP_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

#UP10 PROT
prot_UP10_pathway_gonads_1VSALL_all <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP10_list_gonads_1VSALL)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP10_gonads_1VSALL_all <- prot_UP10_pathway_gonads_1VSALL_all %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP10_pathway_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "06_prot_UP10_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP10_gonads_1VSALL_all, file = file.path(outputdirPATH_gonads, "06b_group_UP10_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

################ MG ################  
#ALL DE PROT
prot_DE_pathway_gonads_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% protDE_gonads_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_ALL_gonads_1VSALL_MG <- prot_DE_pathway_gonads_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DE_pathway_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "04_prot_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_ALL_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads,"04b_group_ALL_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

#DOWN PROT
prot_DOWN_pathway_gonads_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN_list_gonads_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN_gonads_1VSALL_MG <- prot_DOWN_pathway_gonads_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN_pathway_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "05_prot_DOWN_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "05b_group_DOWN_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

# #DOWN10 PROT
prot_DOWN10_pathway_gonads_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% DOWN10_list_gonads_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_DOWN10_gonads_1VSALL_MG <- prot_DOWN10_pathway_gonads_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_DOWN10_pathway_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "05_prot_DOWN10_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_DOWN10_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "05b_group_DOWN10_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

#UP PROT
prot_UP_pathway_gonads_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP_list_gonads_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP_gonads_1VSALL_MG <- prot_UP_pathway_gonads_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP_pathway_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "06_prot_UP_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "06b_group_UP_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

#UP10 PROT
prot_UP10_pathway_gonads_1VSALL_MG <- MG_prot_annot_clean %>%
  filter(MG_prot_annot_clean$ID %in% UP10_list_gonads_1VSALL_MG)

#Nombre de fois ou une voie est vue/assignée dans les protéines DE
group_pathway_UP10_gonads_1VSALL_MG <- prot_UP10_pathway_gonads_1VSALL_MG %>%
  group_by(pathway_gamfocyc) %>%
  dplyr::summarise(nb_prot = n()) %>%
  dplyr::arrange(desc(nb_prot))

write.xlsx(x = prot_UP10_pathway_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "06_prot_UP10_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)
write.xlsx(x = group_pathway_UP10_gonads_1VSALL_MG, file = file.path(outputdirPATH_gonads, "06b_group_UP10_DE_pathway_gonads_1VSALL.xlsx"), overwrite=TRUE)

#-------------------------------------------------------------------
#  05.  SAVING DATA
#-------------------------------------------------------------------
save(#sc, n_prot,
     resume_gills_1VS1, 
     protDE_gills_1VS1, DOWN_list_gills_1VS1, UP_list_gills_1VS1, prot_UP_pathway_gills_1VS1_all, prot_DOWN_pathway_gills_1VS1_all, prot_DE_pathway_gills_1VS1_all,
     protDE_gills_1VS1_MG, DOWN_list_gills_1VS1_MG, UP_list_gills_1VS1_MG, prot_UP_pathway_gills_1VS1_MG, prot_DOWN_pathway_gills_1VS1_MG, prot_DE_pathway_gills_1VS1_MG,
     resume_gills_1VSALL, 
     protDE_gills_1VSALL, DOWN_list_gills_1VSALL, UP_list_gills_1VSALL, prot_UP_pathway_gills_1VSALL_all, prot_DOWN_pathway_gills_1VSALL_all, prot_DE_pathway_gills_1VSALL_all, DOWN10_list_gills_1VSALL, UP10_list_gills_1VSALL, prot_UP10_pathway_gills_1VSALL_all, prot_DOWN10_pathway_gills_1VSALL_all,
     protDE_gills_1VSALL_MG, DOWN_list_gills_1VSALL_MG, UP_list_gills_1VSALL_MG, prot_UP_pathway_gills_1VSALL_MG, prot_DOWN_pathway_gills_1VSALL_MG, prot_DE_pathway_gills_1VSALL_MG, DOWN10_list_gills_1VSALL_MG, UP10_list_gills_1VSALL_MG, prot_UP10_pathway_gills_1VSALL_MG, prot_DOWN10_pathway_gills_1VSALL_MG,
     
     resume_caeca_1VS1, 
     protDE_caeca_1VS1, DOWN_list_caeca_1VS1, UP_list_caeca_1VS1, prot_UP_pathway_caeca_1VS1_all, prot_DOWN_pathway_caeca_1VS1_all, prot_DE_pathway_caeca_1VS1_all, 
     protDE_caeca_1VS1_MG, DOWN_list_caeca_1VS1_MG, UP_list_caeca_1VS1_MG, prot_UP_pathway_caeca_1VS1_MG, prot_DOWN_pathway_caeca_1VS1_MG, prot_DE_pathway_caeca_1VS1_MG,
     resume_caeca_1VSALL, 
     protDE_caeca_1VSALL, DOWN_list_caeca_1VSALL, UP_list_caeca_1VSALL, prot_UP_pathway_caeca_1VSALL_all, prot_DOWN_pathway_caeca_1VSALL_all, prot_DE_pathway_caeca_1VSALL_all, DOWN10_list_caeca_1VSALL, UP10_list_caeca_1VSALL, prot_UP10_pathway_caeca_1VSALL_all, prot_DOWN10_pathway_caeca_1VSALL_all,
     protDE_caeca_1VSALL_MG, DOWN_list_caeca_1VSALL_MG, UP_list_caeca_1VSALL_MG, prot_UP_pathway_caeca_1VSALL_MG, prot_DOWN_pathway_caeca_1VSALL_MG, prot_DE_pathway_caeca_1VSALL_MG, DOWN10_list_caeca_1VSALL_MG, UP10_list_caeca_1VSALL_MG, prot_UP10_pathway_caeca_1VSALL_MG, prot_DOWN10_pathway_caeca_1VSALL_MG,
     
     resume_gonads_1VS1, 
     protDE_gonads_1VS1, DOWN_list_gonads_1VS1, UP_list_gonads_1VS1, prot_UP_pathway_gonads_1VS1_all, prot_DOWN_pathway_gonads_1VS1_all, prot_DE_pathway_gonads_1VS1_all, 
     protDE_gonads_1VS1_MG, DOWN_list_gonads_1VS1_MG, UP_list_gonads_1VS1_MG, prot_UP_pathway_gonads_1VS1_MG, prot_DOWN_pathway_gonads_1VS1_MG, prot_DE_pathway_gonads_1VS1_MG,
     resume_gonads_1VSALL, 
     protDE_gonads_1VSALL, DOWN_list_gonads_1VSALL, UP_list_gonads_1VSALL, prot_UP_pathway_gonads_1VSALL_all, prot_DOWN_pathway_gonads_1VSALL_all, prot_DE_pathway_gonads_1VSALL_all, DOWN10_list_gonads_1VSALL, UP10_list_gonads_1VSALL, prot_UP10_pathway_gonads_1VSALL_all, prot_DOWN10_pathway_gonads_1VSALL_all,
     protDE_gonads_1VSALL_MG, DOWN_list_gonads_1VSALL_MG, UP_list_gonads_1VSALL_MG, prot_UP_pathway_gonads_1VSALL_MG, prot_DOWN_pathway_gonads_1VSALL_MG, prot_DE_pathway_gonads_1VSALL_MG, DOWN10_list_gonads_1VSALL_MG, UP10_list_gonads_1VSALL_MG, prot_UP10_pathway_gonads_1VSALL_MG, prot_DOWN10_pathway_gonads_1VSALL_MG,
     file=file.path(outputdirDE, "list_DE_prot.Rdata"))

save(group_pathway_DOWN_gonads_1VS1_MG,
     group_pathway_UP_gonads_1VS1_MG,
     group_pathway_DOWN_gills_1VSALL_all,
     group_pathway_DOWN_caeca_1VSALL_all,
     file=file.path(outputdirDE, "list_DE_pathway_for_radarplots.Rdata"))
