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

##############################################################################
#               I. DATA IMPORTATION                                          #
##############################################################################
#-------------------------------------------------------------------
#  01. DATA IMPORTATION                           
#-------------------------------------------------------------------
## LOAD DATA FROM DE ANALYSIS
# load(file=file.path(outputdirDE, "list_DE_prot.Rdata"), verbose =TRUE)

load(file=file.path(outputdirDE, "list_DE_pathway_for_radarplots.Rdata"), verbose =TRUE)

#male gonads (down)
group_pathway_DOWN_gonads_1VS1_MG
group_pathway_DOWN10_gonads_1VS1_MG <- head(group_pathway_DOWN_gonads_1VS1_MG, n=10) #top10

#female gonads (up)
group_pathway_UP_gonads_1VS1_MG
group_pathway_UP10_gonads_1VS1_MG <- head(group_pathway_UP_gonads_1VS1_MG, n=10) #top10

#gills
group_pathway_DOWN_gills_1VSALL_all
group_pathway_DOWN10_gills_1VSALL_all <- head(group_pathway_DOWN_gills_1VSALL_all, n=10)  #top10

#caeca
group_pathway_DOWN_caeca_1VSALL_all
group_pathway_DOWN10_caeca_1VSALL_all <- head(group_pathway_DOWN_caeca_1VSALL_all, n=10)  #top10

# table correspondances pathway_gamfocyc et superclasses associees (relevees a la main sur gamfocyc)
superclasses <- read_excel(file.path(datadir, "Correspondance_pathway_gamfocyc_to_superclasses.xlsx"))
# superclasses$pathway_gamfocyc <- gsub("/", "-", superclasses$pathway_gamfocyc)
superclasses <- distinct(superclasses) %>%  # supprimer lignes doublons
  dplyr::arrange(superclasses)

print(superclasses[duplicated(superclasses$pathway_gamfocyc), ] ) #afficher pathway_gamfocyc en doublons (si fautes frappes superclasses par exemple)


##############################################################################
#               I. ALL DE PATHWAYS                                           #
##############################################################################
#-------------------------------------------------------------------
#  01. INTERSECTION OF THE 4 ORGANS                           
#-------------------------------------------------------------------
#### ON NE FAIT PAS CAR PAS LA MEME BASE D'EXPRESSION DIFFERENTIELLE

#-------------------------------------------------------------------
#  01. INTERSECTION OF MALE AND FEMALE ORGANS                           
#-------------------------------------------------------------------

###### SIMPLIFIER PATHWAYS EN FILTRANT CEUX QUI ONT LE MEME NB DE PROTEINES DANS MALE ET FEMALE ET CEUX QUI ON MOINS DU DOUBLE DE PROTEINES DANS L'UN OU DANS L'AUTRE
# ON CHERCHE A ILLUSTRER LES VOIES LES PLUS CONTRASTEES

# Croisement des 2 tables pour trouver voies communes
group_pathway_DOWN_gonads_1VS1_MG
group_pathway_UP_gonads_1VS1_MG


mf <- list(
  male_gonads = group_pathway_DOWN_gonads_1VS1_MG$pathway_gamfocyc, 
  female_gonads = group_pathway_UP_gonads_1VS1_MG$pathway_gamfocyc
) # mf = male-female


tiff(filename=file.path(plotdirRadarplot, "02_Radar_plot_male_female_organs.tiff"))
ggvenn(
  mf, 
  fill_color = c("#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

#venn pour avoir liste 
# library(gplots)
mf.table <- venn(mf)
print(mf.table)

# Extraction de la liste des voies
radar_plot_list_mf <- attr(mf.table, "intersections")$`male_gonads:female_gonads`
             
##Radar plot pathway list
radar_plot_list_mf


## liste pathways et nombre proteines pour chaque organe
# testicules
rp_male <- group_pathway_DOWN_gonads_1VS1_MG %>% 
  filter(group_pathway_DOWN_gonads_1VS1_MG$pathway_gamfocyc %in% radar_plot_list_mf ) %>% 
  dplyr::arrange(pathway_gamfocyc) %>% 
  rename(nb_prot = "nb_prot_male_gonads") #ordre alphabetique
# ovaries
rp_female <- group_pathway_UP_gonads_1VS1_MG %>% 
  filter(group_pathway_UP_gonads_1VS1_MG$pathway_gamfocyc %in% radar_plot_list_mf ) %>% 
  dplyr::arrange(pathway_gamfocyc) %>% 
  rename(nb_prot = "nb_prot_female_gonads") #ordre alphabetique



## tableau récapitulatif en mergant les organes
radar_plot_table <- plyr::join_all(list(rp_male,rp_female), by='pathway_gamfocyc', type='left')


## filtrer table selon deux criteres
radar_plot_table_filtree <- radar_plot_table %>% 
  filter(nb_prot_male_gonads-nb_prot_female_gonads!=0) %>% # garder uniquement les voies qui ont un nombre de protéines DE différents
  filter(nb_prot_male_gonads/nb_prot_female_gonads>=2 | nb_prot_male_gonads/nb_prot_female_gonads<= 1/2) # garder uniquement les voies qui ont un rapport de difference de minimum 2 (dans un sens ou dans l'autre)


##################################### TEST POUR DEFINIR TOP 10 OU 20 EN TRIANT TABLE SELON VALEUR DU RAPPORT
# #rapport dans un sens ?
# test <- radar_plot_table %>% 
#   filter(nb_prot_male_gonads-nb_prot_female_gonads!=0) %>% # garder uniquement les voies qui ont un nombre de protéines DE différents
#   filter(nb_prot_male_gonads/nb_prot_female_gonads>=2 | nb_prot_female_gonads/nb_prot_male_gonads>= 2) %>% # garder uniquement les voies qui ont un rapport de difference de minimum 2 (dans un sens ou dans l'autre)
#   mutate(rapport=nb_prot_male_gonads/nb_prot_female_gonads) %>% #sens 1 division
#   dplyr::arrange(desc(rapport))
# 
# #etendue ?
# test3 <- radar_plot_table %>% 
#   filter(nb_prot_male_gonads-nb_prot_female_gonads!=0) %>% # garder uniquement les voies qui ont un nombre de protéines DE différents
#   filter(nb_prot_male_gonads/nb_prot_female_gonads>=2 | nb_prot_female_gonads/nb_prot_male_gonads>= 2) %>% # garder uniquement les voies qui ont un rapport de difference de minimum 2 (dans un sens ou dans l'autre)
#   mutate(val.absolue=abs(nb_prot_female_gonads-nb_prot_male_gonads)) %>% #valeur absolue de l'étendue (valeur 1 - valeur 2)
#   dplyr::arrange(desc(val.absolue))
####################################################################


## vérifier ce qui a été supprimé
anti_join(radar_plot_table, radar_plot_table_filtree)

##ajouter superclasses des pathways avant exportation
radar_plot_table_finale <- left_join(radar_plot_table_filtree,superclasses,by="pathway_gamfocyc") %>%
  dplyr::arrange(superclasses)

## exporter table pour faire radar plot sur excel
write.xlsx(x = radar_plot_table_finale, file = file.path(outputdirRadarplot, "Radar_plot_table_male_and_female_gonads_rapport2.xlsx"), overwrite=TRUE)



#-------------------------------------------------------------------
#  01. GILLS ONLY                           
#-------------------------------------------------------------------

#gills
View(group_pathway_DOWN_gills_1VSALL_all)
group_pathway_DOWN10_gills_1VSALL_all <- head(group_pathway_DOWN_gills_1VSALL_all, n=11)  #top11 #aller jusqu'à 11 car nb_prot ex aequo sur la 10eme voies


##ajouter superclasses des pathways avant exportation
gills_table_finale <- left_join(group_pathway_DOWN10_gills_1VSALL_all,superclasses, by="pathway_gamfocyc") %>%
  dplyr::arrange(superclasses)


write.xlsx(x = gills_table_finale, file = file.path(outputdirRadarplot, "Barplot_gills_top10.xlsx"), overwrite=TRUE)


##BAR PLOT GILLS

# install.packages("forcats")
library(forcats) #package pour la focntion fct_rev() et fct_infreq() pour réordonner les données



# ggplot(data=gills_table_finale, aes(x=fct_rev(fct_infreq(pathway_gamfocyc)), y=pathway_gamfocyc)) + #beug sur ordre bar plot
#   # geom_bar(stat="identity", aes(fill=superclasses), width=0.4) +
#   geom_bar(stat="identity", aes(fill=superclasses)) +
#   # coord_flip() +
#   theme_minimal() +
#   geom_text(aes(label=nb_prot), vjust=1.6, color="white", size=5.8) +
#   scale_fill_brewer(palette="Dark2") +
#   labs(title = "Gills DE pathways", x="",  y="Number of proteins") +
#   theme(axis.text.x = element_text(angle = 90, hjust=1, size = 16),
#         axis.text.y = element_blank(),
#         plot.title = element_text(size=14),
#         axis.title.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         legend.position = c(0.8, 0.7),
#         legend.title = element_text(size=14),
#         legend.text = element_text(size = 11),
#         legend.background = element_rect(fill = "white"),
#         legend.key.size = unit(0.5, "cm"),
#         legend.key.width = unit(0.5, "cm")) #+
#   # scale_x_discrete(limits=)


ggplot(data=gills_table_finale, aes(x=pathway_gamfocyc, y=nb_prot)) +
  # geom_bar(stat="identity", aes(fill=superclasses), width=0.4) +
  geom_bar(stat="identity", aes(fill=superclasses)) +
  # coord_flip() +
  theme_minimal() +
  geom_text(aes(label=nb_prot), vjust=1.6, color="white", size=5.8) +
  scale_fill_brewer(palette="Dark2") +
  labs(title = "Gills DE pathways", x="",  y="Number of proteins") +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 16),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "top",
        legend.title = element_text(size=14),
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill = "white"),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_line(linetype="blank"),
        # panel.grid.minor = element_line(linetype="blank")
        )


ggsave(filename=file.path(plotdirBarplot, "barplot_top10_gills.tiff"), device="tiff", width=12, height=8, dpi=300)

#-------------------------------------------------------------------
#  01. CAECA ONLY                           
#-------------------------------------------------------------------

#caeca
View(group_pathway_DOWN_caeca_1VSALL_all)
group_pathway_DOWN10_caeca_1VSALL_all <- head(group_pathway_DOWN_caeca_1VSALL_all, n=10)  #top10 
#voir si on arête top10 à 11 pour faire comme gills ?

##ajouter superclasses des pathways avant exportation
caeca_table_finale <- left_join(group_pathway_DOWN10_caeca_1VSALL_all,superclasses,by="pathway_gamfocyc") %>%
  dplyr::arrange(superclasses)


write.xlsx(x = caeca_table_finale, file = file.path(outputdirRadarplot, "Barplot_caeca_top10.xlsx"), overwrite=TRUE)


##BAR PLOT CAECA
# install.packages("forecats")
library(forcats)

# ggplot(data=caeca_table_finale, aes(x=fct_rev(fct_infreq(pathway_gamfocyc)), y=pathway_gamfocyc)) + #beug ordre bar plot
#   geom_bar(stat="identity", aes(fill=superclasses)) +
#   # coord_flip() +
#   theme_minimal() +
#   geom_text(aes(label=nb_prot), vjust=1.6, color="white", size=5.8) +
#   scale_fill_brewer(palette="Set1") +
#   labs(title = "Caeca DE pathways", x="", y="Number of proteins") +
#   theme(axis.text.x = element_text(angle = 90, hjust=1, size = 16),
#         axis.text.y = element_blank(),
#         plot.title = element_text(size=14),
#         axis.title.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         legend.position = c(0.8, 0.7),
#         legend.title = element_text(size=14),
#         legend.text = element_text(size = 11),
#         legend.background = element_rect(fill = "white"),
#         legend.key.size = unit(0.5, "cm"),
#         legend.key.width = unit(0.5, "cm"))

ggplot(data=caeca_table_finale, aes(x=pathway_gamfocyc, y=nb_prot)) +
  geom_bar(stat="identity", aes(fill=superclasses)) +
  # coord_flip() +
  theme_minimal() +
  geom_text(aes(label=nb_prot), vjust=1.6, color="white", size=5.8) +
  scale_fill_brewer(palette="Set1") +
  labs(title = "Caeca DE pathways", x="", y="Number of proteins") +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 16),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "top",
        legend.title = element_text(size=14),
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill = "white"),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_line(linetype="blank"),
        # panel.grid.minor = element_line(linetype="blank")
        )

ggsave(filename=file.path(plotdirBarplot, "barplot_top10_caeca.tiff"), device="tiff", width=12, height=8, dpi=300)



###########################
save(group_pathway_DOWN_gonads_1VS1_MG,
     group_pathway_UP_gonads_1VS1_MG,
     radar_plot_table_finale,
     group_pathway_DOWN_gills_1VSALL_all,
     gills_table_finale,
     group_pathway_DOWN_caeca_1VSALL_all,
     caeca_table_finale,
     file=file.path(outputdirRadarplot, "data_for_patways_analysis.Rdata"))
        