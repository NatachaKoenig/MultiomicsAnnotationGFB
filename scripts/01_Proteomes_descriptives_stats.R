#####################################################################
## PROJECT : GamfoCyc                                              ##
## STUDIES : APPROVE                                               ##
## AUTHOR : Amélie Lafont                                          ##
## DATE : July 2022                                                ##
## SCRIPT : Shotgun proteomics MG/ML validation from GamfoCyc      ##
#####################################################################
#-------------------------------------------------------------------
#  INTRODUCTORY NOTE                            
#-------------------------------------------------------------------
# This script allows to retrieve the attributes (sequence, name ...) 
# of the proteins annotated by CycADS (GamfoCyc database) as belonging
# to the Lipid Metabolism (MG)

#-------------------------------------------------------------------
#  INDEX.SUMMARY OF THE PIPELINE                   
#-------------------------------------------------------------------

#      PACKAGES & SETTINGS 
#      DIRECTORIES  

# I. PROTEOMICS
#   01. DATA IMPORTATION
#   02. RECUPERATION DINFORMATIONS GENERALES SUR LA TABLE S2
#   03. NETTOYAGE DES ACCESSION DANS TABLES s2_GFB

# II. TRANSCRIPTOMICS
#   01. DATA IMPORTATION

# SAVING DATA  

#-------------------------------------------------------------------
#  PACKAGES & SETTINGS                          
#-------------------------------------------------------------------
## Required packages
# Installs missing libraries !
list.of.packages <- c("plyr", "dplyr", "ggplot2", "grid", "gridExtra", "RColorBrewer", "Biobase", "stringr", "dplyr", "tibble", "readxl", "tidyr", "scales") #list of packages required
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
library(scales) #formatte les etiquettes en pourcentage

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



##############################################################################
#                               I. PROTEOMICS                                #
##############################################################################

#-------------------------------------------------------------------
#  01. DATA IMPORTATION                           
#-------------------------------------------------------------------
### Import the s2 proteins table 
raw_s2_gfb <- read.delim(file.path(datadir, "SC_full_names.txt"), 
                         check.names=FALSE, stringsAsFactors=FALSE)
raw_s2_gfb <- as_tibble(raw_s2_gfb) #conversion en tibble


# Import the s2 phenotype table
pData <- read.delim(file = file.path(datadir, "pData_gonad_fr.txt"))

#-------------------------------------------------------------------
#  02. NETTOYAGE DES ACCESSION DANS TABLES s2_GFB
#-------------------------------------------------------------------
s2_gfb_counts <- as.data.frame(raw_s2_gfb)

#nettoyage des noms des échantillons
raw_ech_MS_counts <- str_extract_all(colnames(raw_s2_gfb), "E[0-9]{1,5}")
ech_MS_counts <- unlist(raw_ech_MS_counts)

for (i in 2:ncol(raw_s2_gfb)) {
  colnames(s2_gfb_counts)[i] <- ech_MS_counts[i-1]
}

#nettoyage des noms des protéines
raw_id_MS = str_extract_all(raw_s2_gfb$Accession, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}") #extension identifiant de isoforme
id_MS <- unlist(raw_id_MS)

s2_gfb_counts$Accession <- gsub("TRINITY_", "", raw_s2_gfb$Accession)
s2_gfb_prot_list_unique <- s2_gfb_counts$Accession

s2_gfb_counts <- data.frame(s2_gfb_counts, row.names = 1)

# table des infos s2 apres clean
write.table(s2_gfb_counts, file=(file.path(outputdir, "s2_gfb_counts.txt")), sep = "\t", row.names = FALSE, quote = FALSE) #mettre row.names=TRUE sinon pas de rownames dans le fichier de sortie ?

#-------------------------------------------------------------------
#  03. RECUPERATION DINFORMATIONS GENERALES SUR LA TABLE S2                          
#-------------------------------------------------------------------
### Format S2 table for downstream analysis
s2_gfb <- data.frame(matrix(ncol = 2, nrow = 0)) #creation table qui sera remplie par la boucle juste après
colnames(s2_gfb) <-c("Accession", "Origine_Query")
#Boucle qui récupére les echantillons pour lesquels un contig est comptabilisé
for (i in 1:nrow(raw_s2_gfb)) {
  for (j in 2:ncol(raw_s2_gfb)) {
    if (raw_s2_gfb[i,j] > 0) {
      s2_gfb[nrow(s2_gfb)+1,] <- c(raw_s2_gfb[i,1], names(raw_s2_gfb[j]))
    }
  }
}

### TABLE s2
# extraction code echantillon
raw_ech_MS <- str_extract_all(s2_gfb$Origine_Query, "E[0-9]{1,5}")
ech_MS <- unlist(raw_ech_MS)
length(unique(ech_MS)) # 36 echantillons distincts

# extraction organe
raw_organ_MS <- sapply(strsplit(s2_gfb$Origine_Query, "_"), "[", 5)
organ_MS <- sapply(strsplit(raw_organ_MS, "-"), "[", 1)
length(unique(organ_MS)) # 7 organes distincts

# extraction de ID complet
raw_id_MS = str_extract_all(s2_gfb$Accession, "DN[0-9]*_c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}") #extension identifiant de isoforme
id_MS <- unlist(raw_id_MS)

# extraction du DNxxx
raw_DN_MS <- str_extract_all(s2_gfb$Accession, "DN[0-9]*")
DN_MS <- unlist(raw_DN_MS)

# extraction extension identifiant isoforme
raw_iso_MS = str_extract_all(s2_gfb$Accession, "c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}") #extension identifiant de isoforme
iso_MS <- unlist(raw_iso_MS)

# table des infos s2 apres clean
s2_gfb_clean <- data.frame(ID = id_MS, DN = DN_MS, isoforme = iso_MS, Echantillon = ech_MS, Organ = organ_MS)
write.table(s2_gfb_clean, file=(file.path(outputdir, "s2_gfb_clean.txt")), sep = "\t", row.names = FALSE, quote = FALSE)



##############################################################################
#                             II. TRANSCRIPTOMICS                            #
##############################################################################

##################################### MG #####################################

#-------------------------------------------------------------------
#  01. DATA IMPORTATION OF MG PROTEINS AND PATHWAYS OF GFB                          
#-------------------------------------------------------------------
# Import MG proteins from GamfoCyc
raw_gamfo_MG <- read_excel(file.path(datadir, "v1_MG_pathways_GamfoCyc_GFB.xlsx"), col_names = TRUE)
raw_gamfo_MG <- as_tibble(raw_gamfo_MG) #conversion en tibble

# Nombre de proteines uniques de raw_gamfo_MG (suppression redondance)
length(unique(raw_gamfo_MG$ID_prot)) #nombre de proteines uniques : 4033
gamfo_prot_list_unique_MG <- unique(raw_gamfo_MG$ID_prot)

#-------------------------------------------------------------------
#  02. NB DE PROT ASSOCIEES PAR PATHWAY DU MG POUR GFB (REDUNDANT)                          
#-------------------------------------------------------------------
### Table du nombre de proteines (non unique) pour chaque pathway
#Une prot peut donc se retrouver annotée ou non dans plusieurs voies donc etre presentes sur plusieurs ligne de la table, 
#c'est pour ça que la première voie contient beaucoup de proteine, parce que quand la fonction distincte est applique, elle 
#choisi la première occurrence de la voie associée à la protéine
Nb_prot_each_MG_pathway_redundant <- raw_gamfo_MG %>% 
  group_by(Pathway) %>%
  summarise(nb_of_protein = n_distinct(ID_prot))
#export
write.table(Nb_prot_each_MG_pathway_redundant, file=(file.path(outputdir, "01_Nb_prot_each_MG_pathway_redundant_GFB.txt")), sep = "\t", row.names = FALSE, quote = FALSE)


#-------------------------------------------------------------------
#  A. EXTRACTION LISTE DES DNxxxx DANS raw_gamfo_MG                  
#-------------------------------------------------------------------
# Extraction id complet de la proteine: gsub
id_MG <- gsub("-PA", "", raw_gamfo_MG$ID_prot) #suppression "-PA" dans ID_prot

# extraction du DNxxx
raw_DN_MG <- str_extract_all(id_MG, "DN[0-9]*")
DN_MG <- unlist(raw_DN_MG)

# extraction extension identifiant isoforme
raw_iso_MG = str_extract_all(id_MG, "c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}") #extension identifiant de isoforme
iso_MG <- unlist(raw_iso_MG)

# table des proteines du MG apres clean
MG_prot_clean <- data.frame(ID=id_MG, DN = DN_MG, isoforme = iso_MG, pathway_gamfocyc = raw_gamfo_MG$Pathway)

# export table des proteines du MG selon GamfoCyc
write.table(MG_prot_clean, file=(file.path(outputdir, "03_List_MG_prot_nr_GFB_gamfocyc.txt")), sep = "\t", row.names = FALSE, quote = FALSE)


#-------------------------------------------------------------------
#  B. FILTRER PROTEINES s2 SELON ID PROTEINES MG GAMFOCYC                   
#-------------------------------------------------------------------
### Tri strict sur ID
gamfo_prot_list_MG <- MG_prot_clean$ID
length(gamfo_prot_list_MG) #longueur list prot ML gamfocyc : 10801 proteines

s2_gfb_gamfo_MG_clean <- filter(s2_gfb_clean, s2_gfb_clean$ID %in% gamfo_prot_list_MG) 
# nombre de proteines uniques restantes
length(unique(s2_gfb_gamfo_MG_clean$ID)) # 858 ID  sur les xxx prot du départ
filt_prot_gamfo_MG_unique <- unique(s2_gfb_gamfo_MG_clean$ID)

# nombre echantillons unique restants
length(unique(s2_gfb_gamfo_MG_clean$Echantillon)) # 36 échantillons



#-------------------------------------------------------------------
#  C. AJOUT ANNOTATIONS DANS MG_prot_clean                        
#-------------------------------------------------------------------
raw_MG_prot_annot <- read_excel(file.path(datadir, "out.emapper.annotations.xlsx"), col_names = TRUE, skip=2)
raw_MG_prot_annot<- as_tibble(raw_MG_prot_annot)
MG_prot_annot <- dplyr::select(raw_MG_prot_annot,query,Description)
MG_prot_annot$query <- gsub("TRINITY_", "", MG_prot_annot$query)
MG_prot_annot <- MG_prot_annot %>%
  dplyr::rename(ID = query)
MG_prot_annot

MG_prot_annot_clean <- dplyr::left_join(MG_prot_clean, MG_prot_annot, by ="ID")


##################################### ML #####################################

#-------------------------------------------------------------------
#  01. DATA IMPORTATION OF ML PROTEINS AND PATHWAYS OF GFB                          
#-------------------------------------------------------------------
# Import ML proteins from GamfoCyc
raw_gamfo_ML <- read_excel(file.path(datadir, "v1_ML_pathways_GamfoCyc_GFB.xlsx"), col_names = TRUE)
raw_gamfo_ML <- as_tibble(raw_gamfo_ML) #conversion en tibble

# Nombre de proteines uniques de raw_ML_prot (suppression redondance)
length(unique(raw_gamfo_ML$ID_prot)) #nombre de proteines uniques : 866
gamfo_prot_list_unique_ML <- unique(raw_gamfo_ML$ID_prot)


#-------------------------------------------------------------------
#  02. NB DE PROT ASSOCIEES PAR PATHWAY DU ML POUR GFB (REDUNDANT)                          
#-------------------------------------------------------------------
### Table du nombre de proteines (non unique) pour chaque pathway
#Une prot peut donc se retrouver annotée ou non dans plusieurs voies donc etre presentes sur plusieurs ligne de la table, c'est pour ça que la première voie contient beaucoup de proteine, parce que quand la focntion distincte es tapplique, elle choisi la première occurrence de la voie associée à la protéine
Nb_prot_each_ML_pathway_redundant <- raw_gamfo_ML %>% 
  group_by(Pathway) %>%
  summarise(nb_of_protein = n_distinct(ID_prot))
#export
write.table(Nb_prot_each_ML_pathway_redundant, file=(file.path(outputdir, "01_Nb_prot_each_ML_pathway_redundant_GFB.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

#-------------------------------------------------------------------
#  A. EXTRACTION LISTE DES DNxxxx DANS raw_gamfo_ML                  
#-------------------------------------------------------------------
# Extraction id complet de la proteine: gsub
id_ML <- gsub("-PA", "", raw_gamfo_ML$ID_prot) #suppression "-PA" dans ID_prot

# extraction du DNxxx
raw_DN_ML <- str_extract_all(id_ML, "DN[0-9]*")
DN_ML <- unlist(raw_DN_ML)

# extraction extension identifiant isoforme
raw_iso_ML = str_extract_all(id_ML, "c[0-9]{1,3}_g[0-9]{1,3}_i[0-9]{1,3}\\.p[0-9]{1,3}") #extension identifiant de isoforme
iso_ML <- unlist(raw_iso_ML)

# table des proteines du ML apres clean
ML_prot_clean <- data.frame(ID=id_ML, DN = DN_ML, isoforme = iso_ML, pathway_gamfocyc = raw_gamfo_ML$Pathway)

# export table des proteines du ML selon GamfoCyc
write.table(ML_prot_clean, file=(file.path(outputdir, "03_List_ML_prot_nr_GFB_gamfocyc.txt")), sep = "\t", row.names = FALSE, quote = FALSE)


#-------------------------------------------------------------------
#  B. FILTRER PROTEINES s2 SELON ID PROTEINES ML GAMFOCYC                   
#-------------------------------------------------------------------
### Tri strict sur ID
gamfo_prot_list_ML <- ML_prot_clean$ID
length(gamfo_prot_list_ML) #longueur list prot ML gamfocyc : 2899 proteines

s2_gfb_gamfo_ML_clean <- filter(s2_gfb_clean, s2_gfb_clean$ID %in% gamfo_prot_list_ML) 
# nombre de proteines uniques restantes
length(unique(s2_gfb_gamfo_ML_clean$ID)) # 150 ID  sur les xxx prot du départ
filt_prot_gamfo_ML_unique <- unique(s2_gfb_gamfo_ML_clean$ID)

# nombre echantillons unique restants
length(unique(s2_gfb_gamfo_ML_clean$Echantillon)) # 35 échantillons (E00385 éliminé)


#-------------------------------------------------------------------
#  C. AJOUT ANNOTATIONS DANS ML_prot_clean                        
#-------------------------------------------------------------------
raw_ML_prot_annot <- read_excel(file.path(datadir, "out.emapper.annotations.xlsx"), col_names = TRUE, skip=2)
raw_ML_prot_annot<- as_tibble(raw_ML_prot_annot)
ML_prot_annot <- dplyr::select(raw_ML_prot_annot,query,Description)
ML_prot_annot$query <- gsub("TRINITY_", "", ML_prot_annot$query)
ML_prot_annot <- ML_prot_annot %>%
  dplyr::rename(ID = query)
ML_prot_annot

ML_prot_annot_clean <- dplyr::left_join(ML_prot_clean, ML_prot_annot, by ="ID")

#-------------------------------------------------------------------
#  SAVING DATA 
#-------------------------------------------------------------------
save(s2_gfb_counts, pData, s2_gfb_prot_list_unique, Nb_prot_each_MG_pathway_redundant, MG_prot_clean, MG_prot_annot_clean, gamfo_prot_list_unique_MG, filt_prot_gamfo_MG_unique, ML_prot_clean, ML_prot_annot_clean, gamfo_prot_list_unique_ML, filt_prot_gamfo_ML_unique, file=file.path(outputdir,"data.Rdata"))

