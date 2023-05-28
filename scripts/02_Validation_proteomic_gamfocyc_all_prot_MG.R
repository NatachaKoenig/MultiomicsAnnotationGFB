#####################################################################
## PROJECT : GamfoCyc                                              ##
## STUDIES : APPROVE                                               ##
## AUTHOR : Amélie Lafont                                          ##
## DATE : July 2022                                                ##
## SCRIPT : Shotgum proteomics MG/ML validation from GamfoCyc      ##
#####################################################################
#-------------------------------------------------------------------
#  INTRODUCTORY NOTE                            
#-------------------------------------------------------------------
# This script allows to retrieve the attibutes (sequence, name ...) 
# of the proteins annotated by CycADS (GamfoCyc database) as belonging
# to the Lipid Metabolism (MG)

#-------------------------------------------------------------------
#  INDEX.SUMMARY OF THE PIPELINE                   
#-------------------------------------------------------------------

#      PACKAGES & SETTINGS 
#      DIRECTORIES  

# 01.  DATA IMPORTATION
# 02.  STATISTIQUES PEPTIDES - REPARTITIONS  

#-------------------------------------------------------------------
#  PACKAGES & SETTINGS                          
#-------------------------------------------------------------------
## Required packages
# Installs missing libraries !
list.of.packages <- c("plyr", "dplyr", "ggplot2", "grid", "gridExtra", "RColorBrewer", "Biobase", "stringr", "dplyr", "tibble", "readxl", "tidyr", "scales", "colorspace") #list of packages required
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] #list of packages non installed
if(length(new.packages)) install.packages(new.packages, repos='https://cran.rstudio.com/') #install packages if the new.packages list is not empty

# Installation Biobase
#BiocManager::install("Biobase")

#tools
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(Biobase)
library(stringr)
library(dplyr)
library(tibble)
library(readxl)
library(tidyr)
library(scales)
library(colorspace)#formatte les etiquettes en pourcentage

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
outputdir <- file.path(wdir, "output")

#-------------------------------------------------------------------
#  DATA IMPORTATION
#-------------------------------------------------------------------
load(file=file.path(outputdir,"data.Rdata"), verbose =TRUE)

#-------------------------------------------------------------------
#  STATISTIQUES PROTEINES - REPARTITIONS                   
#-------------------------------------------------------------------

# nombre de proteines uniques restantes
length(unique(s2_gfb_gamfo_MG_clean$ID)) # 858 ID  sur les xxx prot du départ

# nombre echantillons unique restants
length(unique(s2_gfb_gamfo_MG_clean$Echantillon)) # 36 échantillons

########------------ PROTEINES REPARTITION

s2_gfb_gamfo_MG_clean$Organ <- gsub("Testicule", "Gonades", s2_gfb_gamfo_MG_clean$Organ) #modification pour la légende de Testicule en Gonades
s2_gfb_gamfo_MG_clean$Organ <- gsub("Ovaire", "Gonades", s2_gfb_gamfo_MG_clean$Organ) #Ovaire devient Gonades

### répartition protéines selon organes
repartition_prot = table(s2_gfb_gamfo_MG_clean$Organ)
repartition_prot

#Branchie    Caecum  Cephalon  Intestin    Ovaire     Reste Testicule 
#   2677       1649      2524       803       847      2296       914

# Branchie   Caecum Cephalon  Gonades Intestin    Reste 
# 2677     1649     2524     1761      803     2296 

## plot répartition protéines selon organes et echantillon GFB
colourCount <- length(unique(s2_gfb_gamfo_MG_clean$Echantillon)) #nombre de couleur nécessaire
# mycolors <- colorRampPalette(brewer.pal(6,"Blue-Red"))(colourCount) #extension de la palette
mycolors <- diverge_hcl(colourCount, palette = "Green-Brown")

repart_gfb_mg <- ggplot(s2_gfb_gamfo_MG_clean) +
  geom_bar(aes(x=Organ, fill= Echantillon), position="dodge") +
  theme_classic()+
  theme(legend.position = "right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.2,"cm"),
        title = element_text(face="plain", size=20),
        axis.text.x=element_text(face="plain", size=20),
        axis.text.y=element_text(face="plain", size=20),
        axis.title.x=element_text(face="plain", size=20),
        axis.title.y=element_text(face="plain", size=20)) +
  ggtitle("Répartition des protéines présentes dans des voies métaboliques de GfB selon organes et échantillons") +
  xlab("Organes") +
  ylab("Nombre de protéines") +
  scale_fill_manual(values=mycolors, name="Echantillons") +
  scale_y_continuous(limits = c(0,600))

repart_gfb_mg

ggsave(repart_gfb_mg, filename=file.path(plotdir, "03_Number_prot_each_organ_ech_gfb_MG_gonads.tiff"), width=18, height=8, dpi=300)

