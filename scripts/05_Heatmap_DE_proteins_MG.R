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
list.of.packages <- c("plyr", "dplyr", "ggplot2", "grid", "gridExtra", "mixOmics", "minfi", "lumi", "stats", "limma", "edgeR", "Heatplus", "made4", "RColorBrewer", "Biobase", "stringr", "readxl", "tibble", "stringr", "tidyr", "pheatmap") #list of packages required
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


#-------------------------------------------------------------------
#  LOAD DATA                          
#-------------------------------------------------------------------
load(file=file.path(outputdirDE, "data_DE.Rdata"), verbose =TRUE)


#-------------------------------------------------------------------
#  pData EN ANGLAIS POUR ARTICLE                 
#------------------------------------------------------------------
pData_eng <- read.delim(file = file.path(datadir, "pData.txt"))


###############################################################################
################################ HEATMAPS ALL PROT yN2 ########################
###############################################################################
#-------------------------------------------------------------------
#  EXTRACTION SC POUR HEATMAP                 
#------------------------------------------------------------------
yN2counts <- as_tibble(yN2$counts,rownames=NA)


#------------------------------------------------------------------
#  ANNOTATION
#------------------------------------------------------------------
head(yN2counts)
dim(yN2counts) #[1] 2190   36

tyN2 = t(yN2counts)

pData2 <- pData_eng %>%
  dplyr::arrange(SexOrgan, .by_group = TRUE)
annot = data.frame("Organs"=pData2$SexOrgan)
rownames(annot) = rownames(pData2)

lab_col = rownames(annot)
lab_row = colnames(tyN2)
ann_colors = list(Organs = c(Gills_F= "#32B866", Gills_M= "#8FE0AE", Caeca_F="#3498DB", Caeca_M="#96CAED", Cephalon_F="#A624DD", Cephalon_M="#C573E8", 
                              Ovaries_F="#CB001B", Testes_M="#F6798A", Gut_F="#FF9317", Gut_M="#FFBC70", Rest_F="#3A2865", Rest_M="#543A92" ))
my_col_order <- unlist(row.names(pData2))

#-------------------------------------------------------------------
#  HEATMAP TOUTE LES PROTEINES              
#------------------------------------------------------------------
#heatmap all organs and proteins 

tiff(filename=file.path(plotdirDE, "04_heatmap_all_prot_before_DE.tiff"))
pheatmap(yN2counts[,my_col_order],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale="row", 
         main=paste("Enzymes seen in proteomics in GfB organs, filtered \n for (>3 SC & seen in 3 samples)"),
         display_numbers = FALSE,
         fontsize_row = 6,
         fontsize_col = 8,
         fontsize = 8,
         show_colnames = T,
         show_rownames = F,
         annotation_col = annot,
         annotation_colors = ann_colors,
         annotation_legend=T,
         # labels_col = lab_col,
         # labels_row = lab_row,
         # cellheight=4,
         # cellwidth=30,
         # border_color=NA,
         color = colorRampPalette(c("darkgreen", "white", "darkorange"))(n = 200),
         # cutree_cols=6,
         cutree_rows=5
)
dev.off()



###############################################################################
################################ HEATMAPS DE PROT : yN2 #######################
###############################################################################


################################# BRANCHIES ###################################

#-------------------------------------------------------------------
#  EXTRACTION SC POUR HEATMAP                 
#------------------------------------------------------------------
yN2counts <- as_tibble(yN2$counts,rownames=NA)

data_gills <- yN2counts %>%
  filter(rownames(yN2$counts) %in% rownames(all_prot_DE_gills_1VSALL_MG))

#-------------------------------------------------------------------
#  HEATMAP A PARTIR DE yN2                 
#------------------------------------------------------------------
head(data_gills)
dim(data_gills) #[1] 237  36

tdata_gills = t(data_gills)

annot_HM = data.frame("Organs"=pData2$SexOrgan)
rownames(annot_HM) = rownames(pData2)

labels_col = rownames(annot_HM)
labels_row = colnames(tdata_gills)

my_col_order <- unlist(row.names(pData2))

#Heatmap DE prot
tiff(filename=file.path(plotdirDE, "05_heatmap_gills_MG.tiff"))
pheatmap(data_gills[,my_col_order],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale="row", 
         main=paste("Metabolic enzymes validated by proteomics and differentially expressed \n in gills vs. other GfB organs"),
         display_numbers = FALSE,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 8,
         show_colnames = T,
         show_rownames = F,
         # annotation_row = ,
         annotation_col = annot_HM,
         annotation_colors = ann_colors,
         annotation_legend=T,
         # labels_col = labels_col,
         labels_row = labels_row,
         # cellheight=10,
         # cellwidth=30,
         # border_color=NA,
         color = colorRampPalette(c("darkgreen", "white", "darkorange"))(n = 200),
         # cutree_cols=2,
         cutree_rows=2,
)
dev.off()


################################# CAECA ###################################

#-------------------------------------------------------------------
#  EXTRACTION SC POUR HEATMAP                 
#------------------------------------------------------------------
yN2counts <- as_tibble(yN2$counts,rownames=NA)

data_caeca <- yN2counts %>%
  filter(rownames(yN2$counts) %in% rownames(all_prot_DE_caeca_1VSALL_MG))

#-------------------------------------------------------------------
#  HEATMAP A PARTIR DE yN2                 
#------------------------------------------------------------------
head(data_caeca)
dim(data_caeca) #[1] 228  36

tdata_caeca = t(data_caeca)

annot_HM = data.frame("Organs"=pData2$SexOrgan)
rownames(annot_HM) = rownames(pData2)

labels_col = rownames(annot_HM)
labels_row = colnames(tdata_caeca)

my_col_order <- unlist(row.names(pData2))

#Heatmap DE prot
tiff(filename=file.path(plotdirDE, "05_heatmap_caeca_MG.tiff"))
pheatmap(data_caeca[,my_col_order],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale="row", 
         main=paste("Metabolic enzymes validated by proteomics and differentially expressed \n in caeca vs. other GfB organs"),
         display_numbers = FALSE,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 8,
         show_colnames = T,
         show_rownames = F,
         # annotation_row = ,
         annotation_col = annot_HM,
         annotation_colors = ann_colors,
         annotation_legend=T,
         # labels_col = labels_col,
         labels_row = labels_row,
         # cellheight=10,
         # cellwidth=30,
         # border_color=NA,
         color = colorRampPalette(c("darkgreen", "white", "darkorange"))(n = 200),
         # cutree_cols=2,
         cutree_rows=2,
)
dev.off()


################################# GONADES ###################################
################# 1 VS ALL
#-------------------------------------------------------------------
#  EXTRACTION SC POUR HEATMAP                 
#------------------------------------------------------------------
yN2counts <- as_tibble(yN2$counts,rownames=NA)

data_gonads <- yN2counts %>%
  filter(rownames(yN2$counts) %in% rownames(all_prot_DE_gonads_1VSALL_MG))

#-------------------------------------------------------------------
#  HEATMAP A PARTIR DE yN2                 
#------------------------------------------------------------------
head(data_gonads)
dim(data_gonads) #[1] 123  36

tdata_gonads = t(data_gonads)

annot_HM = data.frame("Organs"=pData2$SexOrgan)
rownames(annot_HM) = rownames(pData2)

labels_col = rownames(annot_HM)
labels_row = colnames(tdata_gonads)

my_col_order <- unlist(row.names(pData2))

#Heatmap DE prot
tiff(filename=file.path(plotdirDE, "05_heatmap_gonads_MG_1VSALL.tiff"))
pheatmap(data_gonads[,my_col_order],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale="row", 
         main=paste("Metabolic enzymes validated by proteomics and differentially expressed \n in gonads vs. other GfB organs"),
         display_numbers = FALSE,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 8,
         show_colnames = T,
         show_rownames = F,
         # annotation_row = ,
         annotation_col = annot_HM,
         annotation_colors = ann_colors,
         annotation_legend=T,
         labels_col = labels_col,
         labels_row = labels_row,
         # cellheight=10,
         # cellwidth=30,
         # border_color=NA,
         color = colorRampPalette(c("darkgreen", "white", "darkorange"))(n = 200),
         # cutree_cols=2,
         cutree_rows=2,
)
dev.off()

################# 1VS1 - GONADES F VS GONADES M

#-------------------------------------------------------------------
#  EXTRACTION SC POUR HEATMAP                 
#------------------------------------------------------------------
yN2counts <- as_tibble(yN2$counts,rownames=NA)

data_gonads <- yN2counts %>%
  filter(rownames(yN2$counts) %in% rownames(all_prot_DE_gonads_1VS1_MG))

#-------------------------------------------------------------------
#  HEATMAP A PARTIR DE yN2                 
#------------------------------------------------------------------
head(data_gonads)
dim(data_gonads) #[1] 149  36

tdata_gonads = t(data_gonads)

annot_HM = data.frame("Organs"=pData2$SexOrgan[19:24])
rownames(annot_HM) = rownames(pData2)[19:24]

labels_col = rownames(annot_HM)
labels_row = colnames(tdata_gonads)

ann_colors2 = list(Organs = c(Gonades_F="#CB001B", Gonades_M="#F6798A"))

my_col_order <- unlist(row.names(pData2)[19:24])

#Heatmap DE prot
par(mar= c(6, 6, 0, 0))
tiff(filename=file.path(plotdirDE, "05_heatmap_gonads_MG_1VS1.tiff"), width=240, height=480)
pheatmap(data_gonads[,my_col_order],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale="row", 
         main=paste("Metabolic enzymes validated by proteomics and differentially \n expressed in GfB female gonads vs. male gonads "),
         display_numbers = FALSE,
         fontsize_row = 4,
         fontsize_col = 12,
         fontsize = 6,
         show_colnames = T,
         show_rownames = T,
         # annotation_row = ,
         annotation_col = annot_HM,
         annotation_colors = ann_colors2,
         annotation_legend=T,
         labels_col = labels_col,
         labels_row = labels_row,
         # cellheight=10,
         # cellwidth=30,
         # border_color=NA,
         color = colorRampPalette(c("darkgreen", "white", "darkorange"))(n = 200),
         # cutree_cols=2,
         cutree_rows=2,
)
dev.off()
