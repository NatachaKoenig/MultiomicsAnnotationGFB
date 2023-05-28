#####################################################################
## PROJECT : GamfoCyc                                              ##
## STUDIES : APPROVE                                               ##
## AUTHOR : Natacha Koenig                                         ##
## DATE : April 2023                                               ##
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
library(ggvenn)

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
plotdirRadarplot <- file.path(plotdir, "Radar_plot")
plotdirBarplot <- file.path(plotdir, "Bar_plot")


outputdir <- file.path(wdir, "output")
outputdirDE <- file.path(outputdir, "DE")
outputdirKEGG <- file.path(outputdir, "Mapping_KEGG")
outputdirRadarplot <- file.path(outputdir, "Radar_plot")
outputdirEnrichment<- file.path(outputdir, "Enrichment")

##############################################################################
#               I. DATA IMPORTATION                                          #
##############################################################################
##################################### MG #####################################

#---------------------------------------------------------------------
#  01. DATA IMPORTATION OF ALL PATHWAYS AND NB OF PROT FROM GAMFOCYC                           
#---------------------------------------------------------------------
##### LOAD DATA FROM DESCRIPTIVE ANALYSIS : MG PROTEINS AND PATHWAYS OF GFB (GamfoCyc)
load(file=file.path(outputdir, "data.Rdata"), verbose =TRUE)

head(Nb_prot_each_MG_pathway_redundant)
MG_all_prot_each_pathway <- Nb_prot_each_MG_pathway_redundant
names(MG_all_prot_each_pathway)[1]<-"pathway_gamfocyc"
names(MG_all_prot_each_pathway)[2]<-"nb_prot_before_DE"


############ OU
### Liste des protéines annotées chez gamfocyc et validées par protéomique (n=858)
# head(filt_prot_gamfo_MG_unique)
# TO DO : extraire voies dans lesquelles ces protéines sont impliqquées : voir si cette table n'ewiste pas déjà ?



#### LOAD PATHWAYS FROM PATHWAYS ANALYSIS : PATHWAYS SELECTED FOR RADAR AND BAR PLOTS (script numero 8)
load(file=file.path(outputdirRadarplot, "data_for_patways_analysis.Rdata"), verbose = TRUE)

head(radar_plot_table_finale) #gonads male vs female
head(gills_table_finale) #gills
head(caeca_table_finale) #caeca
##############################################################################
#               I. DATA IMPORTATION                                          #
##############################################################################
#merge table voies selectionnes pour l'analyse des pathways et la table résumant le nombre de protéines par pathways de gamfocyc

#### GONADES
enrichment_table_gonads <- inner_join(radar_plot_table_finale, MG_all_prot_each_pathway, by="pathway_gamfocyc")
#créer nouvelle colonne avec calcul enrichment
enrichment_table_gonads <- mutate(enrichment_table_gonads, enrichment_male_percent = round((nb_prot_male_gonads/nb_prot_before_DE)*100, 2), enrichment_female_percent = round((nb_prot_female_gonads/nb_prot_before_DE)*100, 2))

write.xlsx(x = enrichment_table_gonads, file = file.path(outputdirEnrichment, "Enrichment_table_gonads.xlsx"), overwrite=TRUE)


#### GILLS
enrichment_table_gills <- inner_join(gills_table_finale, MG_all_prot_each_pathway, by="pathway_gamfocyc")
#créer nouvelle colonne avec calcul enrichment
enrichment_table_gills <- mutate(enrichment_table_gills, enrichment_percent = round((nb_prot/nb_prot_before_DE)*100, 2))

write.xlsx(x = enrichment_table_gills, file = file.path(outputdirEnrichment, "Enrichment_table_gills.xlsx"), overwrite=TRUE)


#### CAECA
enrichment_table_caeca <- inner_join(caeca_table_finale, MG_all_prot_each_pathway, by="pathway_gamfocyc")
#créer nouvelle colonne avec calcul enrichment
enrichment_table_caeca <- mutate(enrichment_table_caeca, enrichment_percent = round((nb_prot/nb_prot_before_DE)*100, 2))

write.xlsx(x = enrichment_table_caeca, file = file.path(outputdirEnrichment, "Enrichment_table_caeca.xlsx"), overwrite=TRUE)






