#####################################################################
## PROJECT : GamfoCyc                                              ##
## STUDIES : APPROVE                                               ##
## AUTHOR : Amélie Lafont                                          ##
## DATE : July 2022                                                ##
## SCRIPT : Differential analysis in MG                            ##
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
list.of.packages <- c("plyr", "dplyr", "ggplot2", "grid", "gridExtra", "mixOmics", "minfi", "lumi", "stats", "limma", "edgeR", "Heatplus", "made4", "RColorBrewer", "Biobase", "stringr", "readxl", "tibble", "stringr", "tidyr", "statmod") #list of packages required
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
library(colorspace)
library(Biobase)
library(stringr)
library(readxl)
library(tibble)
library(stringr)
library(tidyr)
library(lumi)
library(statmod)

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
outputdirDE_gills <- file.path(outputdirDE, "gills")
outputdirDE_caeca <- file.path(outputdirDE, "caeca")
outputdirDE_gonads <- file.path(outputdirDE, "gonads")


#-------------------------------------------------------------------
#  LOAD DATA                          
#-------------------------------------------------------------------
load(file=file.path(outputdir,"data.Rdata"), verbose =TRUE)
load(file=file.path(outputdirACP,"data_ACP.Rdata"), verbose =TRUE)

# #############################################################################
#                 I. ANALYSE DIFFERENTIELLE GFB 1 vs 1 ORGANS                 #
# #############################################################################

### EXPLORATION
plotMDS(yN) #avec nom echantillon
points <- c(3,3,18,18,8,8,15,15,16,16,17,17)
colors <- rep(c("#B28634", "#00A495"), 6)

# par(mar= c(5, 6, 4, 2) + 0.1)
par(mar= c(6, 6, 0, 0))
tiff(file=file.path(plotdirDE, "01_MDSplot_SexOrgan.tiff"), width = 480, height = 480)
plotMDS(yN, col=colors[group], pch=points[group], main=paste("Multi-dimensional scaling plot of expression \nprofiles of samples in GfB's proteome"), cex.main=1.8, cex=2, lwd=2, cex.axis=1.5, cex.lab=1.6)
legend("topright", legend=levels(as.factor(group)), pch=points, cex=1.5, col=colors, ncol=2)
dev.off()

# To perform quasi-likelihood F-tests:
design <- model.matrix(~0+group, data=yN$samples)
colnames(design) <- levels(yN$samples$group)
design


#estimation dispersion
yN <- estimateDisp(yN, robust=TRUE)
yN$common.dispersion
# [1] 0.2842571

tiff(file=file.path(plotdirDE, "02_BCV_SexOrgan_1VS1.tiff"))
# pdf(file=file.path(plotdirDE, "02_BCV.pdf"))
plotBCV(yN) #plot estimation dispersion
dev.off()

#ajustement modèle lineraire genralise a chaque proteine
fit <- glmQLFit(yN, design, robust=TRUE) #quasi likeli-hood (QL) pipeline
head(fit$coefficients)

tiff(file=file.path(plotdirDE, "03_QLDisp_SexOrgan_1VS1.tiff"))
# pdf(file=file.path(plotdirDE,"03_QLDisp.pdf"))
plotQLDisp(fit)
dev.off()



#################################### GILLS ####################################  

# #Test de l'expression differentielle
contrast_gills_1VS1 = makeContrasts(Branchies_F - Branchies_M, levels=design) #definition du contraste pour orienter ce qu'on teste (hypothèse nulle)

qlf_gills_1VS1 <- glmQLFTest(fit, contrast=contrast_gills_1VS1) #test DE
topTags(qlf_gills_1VS1) # liste proteine significativement DE

resume_gills_1VS1 = summary(decideTests(qlf_gills_1VS1)) # resume nombre proteine DE (up, down, notsign) avec seuil 5% FDR
resume_gills_1VS1
# 1*Branchies_F -1*Branchies_M
# Down                      1
# NotSig                 2169
# Up                       20

#Plot prot DE : visualisation resultats proteines DE
tiff(file=file.path(plotdirDE, "04_DE_prot_SexOrgan_gills_1VS1.tiff"))
plotMD(qlf_gills_1VS1, cex=1, bg.cex=0.5, hl.col=c("red", "blue"), cex.axis=1, cex.lab=1.5) #plot LFC againt log counts per million avec DE gene highlith
abline(h=c(-1,1), col="blue") #LFC de 2
dev.off()

# Liste prot DE for heatmap
#keeping the protein list according LFC value or FDR value
all_prot_DE_gills_1VS1 <- topTags(qlf_gills_1VS1, n = Inf, adjust.method = "BH", p = 0.05)$table
DE.up_gills_1VS1 <- all_prot_DE_gills_1VS1[all_prot_DE_gills_1VS1$logFC>0, ] #table up reg prot
DE.down_gills_1VS1 <- all_prot_DE_gills_1VS1[all_prot_DE_gills_1VS1$logFC<0, ] #table down reg prot

all_prot_DE_gills_1VS1_MG <- filter(all_prot_DE_gills_1VS1, all_prot_DE_gills_1VS1$genes %in% filt_prot_gamfo_MG_unique)
DE.up_gills_1VS1_MG <- filter(DE.up_gills_1VS1, DE.up_gills_1VS1$genes %in% filt_prot_gamfo_MG_unique)
DE.down_gills_1VS1_MG <- filter(DE.down_gills_1VS1, DE.down_gills_1VS1$genes %in% filt_prot_gamfo_MG_unique)

all_prot_DE_gills_1VS1_ML <- filter(all_prot_DE_gills_1VS1, all_prot_DE_gills_1VS1$genes %in% filt_prot_gamfo_ML_unique)
DE.up_gills_1VS1_ML <- filter(DE.up_gills_1VS1, DE.up_gills_1VS1$genes %in% filt_prot_gamfo_ML_unique)
DE.down_gills_1VS1_ML <- filter(DE.down_gills_1VS1, DE.down_gills_1VS1$genes %in% filt_prot_gamfo_ML_unique)

#all
write.table(all_prot_DE_gills_1VS1, file=file.path(outputdirDE_gills, "01_ALL_prot_DE_gills_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gills_1VS1_MG, file=file.path(outputdirDE_gills, "01_ALL_prot_DE_gills_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gills_1VS1_ML, file=file.path(outputdirDE_gills, "01_ALL_prot_DE_gills_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

# #up
write.table(DE.up_gills_1VS1, file=file.path(outputdirDE_gills, "01_UP_prot_DE_gills_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gills_1VS1_MG, file=file.path(outputdirDE_gills, "01_UP_prot_DE_gills_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gills_1VS1_ML, file=file.path(outputdirDE_gills, "01_UP_prot_DE_gills_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

# #down
write.table(DE.down_gills_1VS1, file=file.path(outputdirDE_gills, "01_DOWN_prot_DE_gills_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gills_1VS1_MG, file=file.path(outputdirDE_gills, "01_DOWN_prot_DE_gills_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gills_1VS1_ML, file=file.path(outputdirDE_gills, "01_DOWN_prot_DE_gills_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

#################################### CAECA ####################################  

# #Test de l'expression differentielle
contrast_caeca_1VS1 = makeContrasts(Caeca_F - Caeca_M, levels=design) #definition du contraste pour orienter ce qu'on teste (hypothèse nulle)

qlf_caeca_1VS1 <- glmQLFTest(fit, contrast=contrast_caeca_1VS1) #test DE
topTags(qlf_caeca_1VS1) # liste proteine significativement DE

resume_caeca_1VS1 = summary(decideTests(qlf_caeca_1VS1)) # resume nombre proteine DE (up, down, notsign) avec seuil 5% FDR
resume_caeca_1VS1
# 1*Caeca_F -1*Caeca_M
# Down                        2
# NotSig                   2155
# Up                         33

#Plot prot DE : visualisation resultats proteines DE
tiff(file=file.path(plotdirDE, "04_DE_prot_SexOrgan_caeca_1VS1.tiff"))
plotMD(qlf_caeca_1VS1, cex=1, bg.cex=0.5, hl.col=c("red", "blue"), cex.axis=1, cex.lab=1.5) #plot LFC againt log counts per million avec DE gene highlith
abline(h=c(-1,1), col="blue") #LFC de 2
dev.off()

# Liste prot DE for heatmap
#keeping the protein list according LFC value or FDR value
all_prot_DE_caeca_1VS1 <- topTags(qlf_caeca_1VS1, n = Inf, adjust.method = "BH", p = 0.05)$table
DE.up_caeca_1VS1 <- all_prot_DE_caeca_1VS1[all_prot_DE_caeca_1VS1$logFC>0, ] #table up reg prot
DE.down_caeca_1VS1 <- all_prot_DE_caeca_1VS1[all_prot_DE_caeca_1VS1$logFC<0, ] #table down reg prot

all_prot_DE_caeca_1VS1_MG <- filter(all_prot_DE_caeca_1VS1, all_prot_DE_caeca_1VS1$genes %in% filt_prot_gamfo_MG_unique)
DE.up_caeca_1VS1_MG <- filter(DE.up_caeca_1VS1, DE.up_caeca_1VS1$genes %in% filt_prot_gamfo_MG_unique)
DE.down_caeca_1VS1_MG <- filter(DE.down_caeca_1VS1, DE.down_caeca_1VS1$genes %in% filt_prot_gamfo_MG_unique)

all_prot_DE_caeca_1VS1_ML <- filter(all_prot_DE_caeca_1VS1, all_prot_DE_caeca_1VS1$genes %in% filt_prot_gamfo_ML_unique)
DE.up_caeca_1VS1_ML <- filter(DE.up_caeca_1VS1, DE.up_caeca_1VS1$genes %in% filt_prot_gamfo_ML_unique)
DE.down_caeca_1VS1_ML <- filter(DE.down_caeca_1VS1, DE.down_caeca_1VS1$genes %in% filt_prot_gamfo_ML_unique)

#all
write.table(all_prot_DE_caeca_1VS1, file=file.path(outputdirDE_caeca, "01_ALL_prot_DE_caeca_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_caeca_1VS1_MG, file=file.path(outputdirDE_caeca, "01_ALL_prot_DE_caeca_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_caeca_1VS1_ML, file=file.path(outputdirDE_caeca, "01_ALL_prot_DE_caeca_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

# #up
write.table(DE.up_caeca_1VS1, file=file.path(outputdirDE_caeca, "01_UP_prot_DE_caeca_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_caeca_1VS1_MG, file=file.path(outputdirDE_caeca, "01_UP_prot_DE_caeca_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_caeca_1VS1_ML, file=file.path(outputdirDE_caeca, "01_UP_prot_DE_caeca_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

# #down
write.table(DE.down_caeca_1VS1, file=file.path(outputdirDE_caeca, "01_DOWN_prot_DE_caeca_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_caeca_1VS1_MG, file=file.path(outputdirDE_caeca, "01_DOWN_prot_DE_caeca_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_caeca_1VS1_ML, file=file.path(outputdirDE_caeca, "01_DOWN_prot_DE_caeca_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

#################################### GONADS ####################################  

# #Test de l'expression differentielle
contrast_gonads_1VS1 = makeContrasts(Gonades_F - Gonades_M, levels=design) #definition du contraste pour orienter ce qu'on teste (hypothèse nulle)

qlf_gonads_1VS1 <- glmQLFTest(fit, contrast=contrast_gonads_1VS1) #test DE
topTags(qlf_gonads_1VS1) # liste proteine significativement DE

resume_gonads_1VS1 = summary(decideTests(qlf_gonads_1VS1)) # resume nombre proteine DE (up, down, notsign) avec seuil 5% FDR
resume_gonads_1VS1
# 1*Gonades_F -1*Gonades_M
# Down                        343
# NotSig                     1660
# Up                          187

#Plot prot DE : visualisation resultats proteines DE
tiff(file=file.path(plotdirDE, "04_DE_prot_SexOrgan_gonads_1VS1.tiff"))
plotMD(qlf_gonads_1VS1, cex=1, bg.cex=0.5, hl.col=c("red", "blue"), cex.axis=1, cex.lab=1.5) #plot LFC againt log counts per million avec DE gene highlith
abline(h=c(-1,1), col="blue") #LFC de 2
dev.off()

# Liste prot DE for heatmap
#keeping the protein list according LFC value or FDR value
all_prot_DE_gonads_1VS1 <- topTags(qlf_gonads_1VS1, n = Inf, adjust.method = "BH", p = 0.05)$table
DE.up_gonads_1VS1 <- all_prot_DE_gonads_1VS1[all_prot_DE_gonads_1VS1$logFC>0, ] #table up reg prot
DE.down_gonads_1VS1 <- all_prot_DE_gonads_1VS1[all_prot_DE_gonads_1VS1$logFC<0, ] #table down reg prot

all_prot_DE_gonads_1VS1_MG <- filter(all_prot_DE_gonads_1VS1, all_prot_DE_gonads_1VS1$genes %in% filt_prot_gamfo_MG_unique)
DE.up_gonads_1VS1_MG <- filter(DE.up_gonads_1VS1, DE.up_gonads_1VS1$genes %in% filt_prot_gamfo_MG_unique)
DE.down_gonads_1VS1_MG <- filter(DE.down_gonads_1VS1, DE.down_gonads_1VS1$genes %in% filt_prot_gamfo_MG_unique)

all_prot_DE_gonads_1VS1_ML <- filter(all_prot_DE_gonads_1VS1, all_prot_DE_gonads_1VS1$genes %in% filt_prot_gamfo_ML_unique)
DE.up_gonads_1VS1_ML <- filter(DE.up_gonads_1VS1, DE.up_gonads_1VS1$genes %in% filt_prot_gamfo_ML_unique)
DE.down_gonads_1VS1_ML <- filter(DE.down_gonads_1VS1, DE.down_gonads_1VS1$genes %in% filt_prot_gamfo_ML_unique)

#all
write.table(all_prot_DE_gonads_1VS1, file=file.path(outputdirDE_gonads, "01_ALL_prot_DE_gonads_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gonads_1VS1_MG, file=file.path(outputdirDE_gonads, "01_ALL_prot_DE_gonads_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gonads_1VS1_ML, file=file.path(outputdirDE_gonads, "01_ALL_prot_DE_gonads_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

# #up
write.table(DE.up_gonads_1VS1, file=file.path(outputdirDE_gonads, "01_UP_prot_DE_gonads_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gonads_1VS1_MG, file=file.path(outputdirDE_gonads, "01_UP_prot_DE_gonads_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gonads_1VS1_ML, file=file.path(outputdirDE_gonads, "01_UP_prot_DE_gonads_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

# #down
write.table(DE.down_gonads_1VS1, file=file.path(outputdirDE_gonads, "01_DOWN_prot_DE_gonads_1VS1.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gonads_1VS1_MG, file=file.path(outputdirDE_gonads, "01_DOWN_prot_DE_gonads_1VS1_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gonads_1VS1_ML, file=file.path(outputdirDE_gonads, "01_DOWN_prot_DE_gonads_1VS1_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

###############################################################################
#                II. ANALYSE DIFFERENTIELLE GFB 1 vs ALL ORGANS               #
###############################################################################

group = as.factor(pData$Organ)

yF<-DGEList(counts=s2_gfb_countsF,
            genes=row.names(s2_gfb_countsF),
            group= group)
dim(yF)
head(yF$samples) 
yN <- calcNormFactors(yF) #normalization of the filtered counts 
dim(yN)
head(yN$samples) #containing sample informations

### EXPLORATION DONNEES
# plotMDS(yN) #avec nom echantillon
points <- c(3,18,8,15,16,17)
colors <- palette.colors(6, "Okabe-Ito") #colorblind safe palette

par(mar= c(6, 6, 0, 0))
tiff(file=file.path(plotdirDE, "01_MDSplot_Organ.tiff"), width = 480, height = 480)
plotMDS(yN, col=colors[group], pch=points[group], main=paste("Multi-dimensional scaling to plot expression \nprofiles of samples in GfB's proteome"), cex.main=1.8, cex=2, lwd=2, cex.axis=1.5, cex.lab=2)
legend("topright", legend=levels(as.factor(group)), pch=points, cex=1.5, col=colors, ncol=2)
dev.off()

# To perform quasi-likelihood F-test, we need to define the design matrix (elle vient du design experimental)
design <- model.matrix(~0+group, data=yN$samples)
colnames(design) <- levels(yN$samples$group)
design

#estimation dispersion
yN <- estimateDisp(yN, robust=TRUE)
yN$common.dispersion
# [1] 0.6477869

tiff(file=file.path(plotdirDE, "02_BCV_Organ_1VSALL.tiff"))
# pdf(file=file.path(plotdirDE, "02_BCV.pdf"))
plotBCV(yN) #plot estimation dispersion
dev.off()


#Ajustement du modèle lineaire généralisé (GLM) aux données et au design
fit <- glmQLFit(yN, design, robust=TRUE)
head(fit$coefficients)


tiff(file=file.path(plotdirDE, "03_QLDisp_Organ_1VSALL.tiff"))
# pdf(file=file.path(plotdirDE,"03_QLDisp.pdf"))
plotQLDisp(fit)
dev.off()

################################# GILLS ###################################

# Make Contrast : Control VS The Rest method (dépend de la question que l'on se pose / hypothèse nulle)
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#control-versus-the-rest
contrast_gills_1VSALL = makeContrasts((Cephalon+Caeca+Intestin+Reste+Gonades)/5-(Branchies),levels = colnames(design)) #defini au moment des paramètres
# contrast = makeContrasts((Cephalon_M+Caeca_M+Gut_M+Gonads_M+Rest_M+Cephalon_F+Caeca_F+Gut_F+Gonads_F+Rest_F)/10-(Branchies_M+Branchies_F)/2,levels = colnames(design)) #defini au moment des paramètres

formule_gills_1VSALL = contrast_gills_1VSALL#défini au début du script

#Test de l'expression différentielle
qlf_gills_1VSALL <- glmQLFTest(fit, contrast=makeContrasts(formule_gills_1VSALL, levels = colnames(design))) #Branchie vs Caecum
topTags(qlf_gills_1VSALL)

#Nombre de proteines differentiellement exprimées (DE)
resume_gills_1VSALL = summary(decideTests(qlf_gills_1VSALL)) #total number od DE prot with 5% FDR
resume_gills_1VSALL
# -1*Branchies 0.2*Caeca 0.2*Cephalon 0.2*Gonades 0.2*Intestin 0.2*Reste
# Down                                                                      653
# NotSig                                                                   1397
# Up                                                                        140

#Plot prot DE
tiff(file=file.path(plotdirDE,"04_DE_prot_Organ_gills_1VSALL.tiff"))
plotMD(qlf_gills_1VSALL, cex=1, bg.cex=0.5, hl.col=c("red", "blue"), cex.axis=1, cex.lab=1.5) #plot LFC againt log counts per million avec DE gene highlith
abline(h=c(-1,1), col="blue") #LFC de 2
dev.off()

# Liste prot DE pour la heatmap
#keeping the protein list according LFC value or FDR value 
all_prot_DE_gills_1VSALL <- topTags(qlf_gills_1VSALL, n = Inf, adjust.method = "BH", p = 0.05)$table
DE.up_gills_1VSALL <- all_prot_DE_gills_1VSALL[all_prot_DE_gills_1VSALL$logFC>0, ] #table up reg prot
DE.up.10_gills_1VSALL <- head(all_prot_DE_gills_1VSALL[all_prot_DE_gills_1VSALL$logFC>0, ], n=10) #table up reg prot top 10
DE.down_gills_1VSALL <- all_prot_DE_gills_1VSALL[all_prot_DE_gills_1VSALL$logFC<0, ] #table down reg prot
DE.down.10_gills_1VSALL <- head(all_prot_DE_gills_1VSALL[all_prot_DE_gills_1VSALL$logFC<0, ], n=10) #table down reg prot top10

all_prot_DE_gills_1VSALL_MG <- filter(all_prot_DE_gills_1VSALL, all_prot_DE_gills_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.up_gills_1VSALL_MG <- filter(DE.up_gills_1VSALL, DE.up_gills_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.up.10_gills_1VSALL_MG <- head(all_prot_DE_gills_1VSALL_MG[all_prot_DE_gills_1VSALL_MG$logFC>0, ], n=10) #table up reg prot top 10
DE.down_gills_1VSALL_MG <- filter(DE.down_gills_1VSALL, DE.down_gills_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.down.10_gills_1VSALL_MG <- head(all_prot_DE_gills_1VSALL_MG[all_prot_DE_gills_1VSALL_MG$logFC<0, ], n=10) #table up reg prot top 10

all_prot_DE_gills_1VSALL_ML <- filter(all_prot_DE_gills_1VSALL, all_prot_DE_gills_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.up_gills_1VSALL_ML <- filter(DE.up_gills_1VSALL, DE.up_gills_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.up.10_gills_1VSALL_ML <- head(all_prot_DE_gills_1VSALL_ML[all_prot_DE_gills_1VSALL_ML$logFC>0, ], n=10) #table up reg prot top 10
DE.down_gills_1VSALL_ML <- filter(DE.down_gills_1VSALL, DE.down_gills_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.down.10_gills_1VSALL_ML <- head(all_prot_DE_gills_1VSALL_ML[all_prot_DE_gills_1VSALL_ML$logFC<0, ], n=10) #table up reg prot top 10



#export de ces listes de protéines
#all
write.table(all_prot_DE_gills_1VSALL, file=file.path(outputdirDE_gills,"01_ALL_prot_DE_gills_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gills_1VSALL_MG, file=file.path(outputdirDE_gills, "01_ALL_prot_DE_gills_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gills_1VSALL_ML, file=file.path(outputdirDE_gills, "01_ALL_prot_DE_gills_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

#up
write.table(DE.up_gills_1VSALL, file=file.path(outputdirDE_gills,"02_UP_prot_DE_gills_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_gills_1VSALL, file=file.path(outputdirDE_gills,"02_UP_prot_DE_TOP10_gills_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gills_1VSALL_MG, file=file.path(outputdirDE_gills, "02_UP_prot_DE_gills_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_gills_1VSALL_MG, file=file.path(outputdirDE_gills, "02_UP_prot_DE_TOP10_gills_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gills_1VSALL_ML, file=file.path(outputdirDE_gills, "02_UP_prot_DE_gills_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_gills_1VSALL_ML, file=file.path(outputdirDE_gills, "02_UP_prot_DE_TOP10_gills_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


#down
write.table(DE.down_gills_1VSALL, file=file.path(outputdirDE_gills,"02_DOWN_prot_DE_gills_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_gills_1VSALL, file=file.path(outputdirDE_gills,"02_DOWN_prot_DE_TOP10_gills_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gills_1VSALL_MG, file=file.path(outputdirDE_gills, "02_DOWN_prot_DE_gills_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_gills_1VSALL_MG, file=file.path(outputdirDE_gills, "02_DOWN_prot_DE_TOP10_gills_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gills_1VSALL_ML, file=file.path(outputdirDE_gills, "02_DOWN_prot_DE_gills_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_gills_1VSALL_ML, file=file.path(outputdirDE_gills, "02_DOWN_prot_DE_TOP10_gills_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


#################################### CAECA ####################################

# Make Contrast : Control VS The Rest method (dépend de la question que l'on se pose / hypothèse nulle)
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#control-versus-the-rest
contrast_caeca_1VSALL = makeContrasts((Cephalon+Branchies+Intestin+Reste+Gonades)/5-(Caeca),levels = colnames(design)) #defini au moment des paramètres

formule_caeca_1VSALL = contrast_caeca_1VSALL #défini au début du script

#Test de l'expression différentielle
qlf_caeca_1VSALL <- glmQLFTest(fit, contrast=makeContrasts(formule_caeca_1VSALL, levels = colnames(design))) #Branchie vs Caecum
topTags(qlf_caeca_1VSALL)

#Nombre de proteines differentiellement exprimées (DE)
resume_caeca_1VSALL = summary(decideTests(qlf_caeca_1VSALL)) #total number od DE prot with 5% FDR
resume_caeca_1VSALL
# 0.2*Branchies -1*Caeca 0.2*Cephalon 0.2*Gonades 0.2*Intestin 0.2*Reste
# Down                                                                      635
# NotSig                                                                   1265
# Up                                                                        290

#Plot prot DE
tiff(file=file.path(plotdirDE,"04_DE_prot_Organ_caeca_1VSALL.tiff"))
plotMD(qlf_caeca_1VSALL, cex=1, bg.cex=0.5, hl.col=c("red", "blue"), cex.axis=1, cex.lab=1.5) #plot LFC againt log counts per million avec DE gene highlith
abline(h=c(-1,1), col="blue") #LFC de 2
dev.off()

# Liste prot DE pour la heatmap
#keeping the protein list according LFC value or FDR value 
all_prot_DE_caeca_1VSALL <- topTags(qlf_caeca_1VSALL, n = Inf, adjust.method = "BH", p = 0.05)$table
DE.up_caeca_1VSALL <- all_prot_DE_caeca_1VSALL[all_prot_DE_caeca_1VSALL$logFC>0, ] #table up reg prot
DE.up.10_caeca_1VSALL <- head(all_prot_DE_caeca_1VSALL[all_prot_DE_caeca_1VSALL$logFC>0, ], n=10) #table up reg prot top 10
DE.down_caeca_1VSALL <- all_prot_DE_caeca_1VSALL[all_prot_DE_caeca_1VSALL$logFC<0, ] #table down reg prot
DE.down.10_caeca_1VSALL <- head(all_prot_DE_caeca_1VSALL[all_prot_DE_caeca_1VSALL$logFC<0, ], n=10) #table down reg prot top10

all_prot_DE_caeca_1VSALL_MG <- filter(all_prot_DE_caeca_1VSALL, all_prot_DE_caeca_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.up_caeca_1VSALL_MG <- filter(DE.up_caeca_1VSALL, DE.up_caeca_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.up.10_caeca_1VSALL_MG <- head(all_prot_DE_caeca_1VSALL_MG[all_prot_DE_caeca_1VSALL_MG$logFC>0, ], n=10) #table up reg prot top 10
DE.down_caeca_1VSALL_MG <- filter(DE.down_caeca_1VSALL, DE.down_caeca_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.down.10_caeca_1VSALL_MG <- head(all_prot_DE_caeca_1VSALL_MG[all_prot_DE_caeca_1VSALL_MG$logFC<0, ], n=10) #table up reg prot top 10

all_prot_DE_caeca_1VSALL_ML <- filter(all_prot_DE_caeca_1VSALL, all_prot_DE_caeca_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.up_caeca_1VSALL_ML <- filter(DE.up_caeca_1VSALL, DE.up_caeca_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.up.10_caeca_1VSALL_ML <- head(all_prot_DE_caeca_1VSALL_ML[all_prot_DE_caeca_1VSALL_ML$logFC>0, ], n=10) #table up reg prot top 10
DE.down_caeca_1VSALL_ML <- filter(DE.down_caeca_1VSALL, DE.down_caeca_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.down.10_caeca_1VSALL_ML <- head(all_prot_DE_caeca_1VSALL_ML[all_prot_DE_caeca_1VSALL_ML$logFC<0, ], n=10) #table up reg prot top 10



#export de ces listes de protéines
#all
write.table(all_prot_DE_caeca_1VSALL, file=file.path(outputdirDE_caeca,"01_ALL_prot_DE_caeca_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_caeca_1VSALL_MG, file=file.path(outputdirDE_caeca, "01_ALL_prot_DE_caeca_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_caeca_1VSALL_ML, file=file.path(outputdirDE_caeca, "01_ALL_prot_DE_caeca_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

#up
write.table(DE.up_caeca_1VSALL, file=file.path(outputdirDE_caeca,"02_UP_prot_DE_caeca_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_caeca_1VSALL, file=file.path(outputdirDE_caeca,"02_UP_prot_DE_TOP10_caeca_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_caeca_1VSALL_MG, file=file.path(outputdirDE_caeca, "02_UP_prot_DE_caeca_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_caeca_1VSALL_MG, file=file.path(outputdirDE_caeca, "02_UP_prot_DE_TOP10_caeca_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_caeca_1VSALL_ML, file=file.path(outputdirDE_caeca, "02_UP_prot_DE_caeca_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_caeca_1VSALL_ML, file=file.path(outputdirDE_caeca, "02_UP_prot_DE_TOP10_caeca_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


#down
write.table(DE.down_caeca_1VSALL, file=file.path(outputdirDE_caeca,"02_DOWN_prot_DE_caeca_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_caeca_1VSALL, file=file.path(outputdirDE_caeca,"02_DOWN_prot_DE_TOP10_caeca_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_caeca_1VSALL_MG, file=file.path(outputdirDE_caeca, "02_DOWN_prot_DE_caeca_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_caeca_1VSALL_MG, file=file.path(outputdirDE_caeca, "02_DOWN_prot_DE_TOP10_caeca_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_caeca_1VSALL_ML, file=file.path(outputdirDE_caeca, "02_DOWN_prot_DE_caeca_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_caeca_1VSALL_ML, file=file.path(outputdirDE_caeca, "02_DOWN_prot_DE_TOP10_caeca_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


#################################### GONADS ####################################

# Make Contrast : Control VS The Rest method (dépend de la question que l'on se pose / hypothèse nulle)
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#control-versus-the-rest
contrast_gonads_1VSALL = makeContrasts((Cephalon+Branchies+Intestin+Reste+Caeca)/5-(Gonades),levels = colnames(design)) #defini au moment des paramètres

formule_gonads_1VSALL = contrast_gonads_1VSALL #défini au début du script

#Test de l'expression différentielle
qlf_gonads_1VSALL <- glmQLFTest(fit, contrast=makeContrasts(formule_gonads_1VSALL, levels = colnames(design))) #Branchie vs Caecum
topTags(qlf_gonads_1VSALL)

#Nombre de proteines differentiellement exprimées (DE)
resume_gonads_1VSALL = summary(decideTests(qlf_gonads_1VSALL)) #total number od DE prot with 5% FDR
resume_gonads_1VSALL
# 0.2*Branchies 0.2*Caeca 0.2*Cephalon -1*Gonades 0.2*Intestin 0.2*Reste
# Down                                                                      413
# NotSig                                                                   1634
# Up                                                                        143

#Plot prot DE
tiff(file=file.path(plotdirDE,"04_DE_prot_Organ_gonads_1VSALL.tiff"))
plotMD(qlf_gonads_1VSALL, cex=1, bg.cex=0.5, hl.col=c("red", "blue"), cex.axis=1, cex.lab=1.5) #plot LFC againt log counts per million avec DE gene highlith
abline(h=c(-1,1), col="blue") #LFC de 2
dev.off()

# Liste prot DE pour la heatmap
#keeping the protein list according LFC value or FDR value 
all_prot_DE_gonads_1VSALL <- topTags(qlf_gonads_1VSALL, n = Inf, adjust.method = "BH", p = 0.05)$table
DE.up_gonads_1VSALL <- all_prot_DE_gonads_1VSALL[all_prot_DE_gonads_1VSALL$logFC>0, ] #table up reg prot
DE.up.10_gonads_1VSALL <- head(all_prot_DE_gonads_1VSALL[all_prot_DE_gonads_1VSALL$logFC>0, ], n=10) #table up reg prot top 10
DE.down_gonads_1VSALL <- all_prot_DE_gonads_1VSALL[all_prot_DE_gonads_1VSALL$logFC<0, ] #table down reg prot
DE.down.10_gonads_1VSALL <- head(all_prot_DE_gonads_1VSALL[all_prot_DE_gonads_1VSALL$logFC<0, ], n=10) #table down reg prot top10

all_prot_DE_gonads_1VSALL_MG <- filter(all_prot_DE_gonads_1VSALL, all_prot_DE_gonads_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.up_gonads_1VSALL_MG <- filter(DE.up_gonads_1VSALL, DE.up_gonads_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.up.10_gonads_1VSALL_MG <- head(all_prot_DE_gonads_1VSALL_MG[all_prot_DE_gonads_1VSALL_MG$logFC>0, ], n=10) #table up reg prot top 10
DE.down_gonads_1VSALL_MG <- filter(DE.down_gonads_1VSALL, DE.down_gonads_1VSALL$genes %in% filt_prot_gamfo_MG_unique)
DE.down.10_gonads_1VSALL_MG <- head(all_prot_DE_gonads_1VSALL_MG[all_prot_DE_gonads_1VSALL_MG$logFC<0, ], n=10) #table up reg prot top 10

all_prot_DE_gonads_1VSALL_ML <- filter(all_prot_DE_gonads_1VSALL, all_prot_DE_gonads_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.up_gonads_1VSALL_ML <- filter(DE.up_gonads_1VSALL, DE.up_gonads_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.up.10_gonads_1VSALL_ML <- head(all_prot_DE_gonads_1VSALL_ML[all_prot_DE_gonads_1VSALL_ML$logFC>0, ], n=10) #table up reg prot top 10
DE.down_gonads_1VSALL_ML <- filter(DE.down_gonads_1VSALL, DE.down_gonads_1VSALL$genes %in% filt_prot_gamfo_ML_unique)
DE.down.10_gonads_1VSALL_ML <- head(all_prot_DE_gonads_1VSALL_ML[all_prot_DE_gonads_1VSALL_ML$logFC<0, ], n=10) #table up reg prot top 10



#export de ces listes de protéines
#all
write.table(all_prot_DE_gonads_1VSALL, file=file.path(outputdirDE_gonads,"01_ALL_prot_DE_gonads_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gonads_1VSALL_MG, file=file.path(outputdirDE_gonads, "01_ALL_prot_DE_gonads_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(all_prot_DE_gonads_1VSALL_ML, file=file.path(outputdirDE_gonads, "01_ALL_prot_DE_gonads_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

#up
write.table(DE.up_gonads_1VSALL, file=file.path(outputdirDE_gonads,"02_UP_prot_DE_gonads_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_gonads_1VSALL, file=file.path(outputdirDE_gonads,"02_UP_prot_DE_TOP10_gonads_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gonads_1VSALL_MG, file=file.path(outputdirDE_gonads, "02_UP_prot_DE_gonads_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_gonads_1VSALL_MG, file=file.path(outputdirDE_gonads, "02_UP_prot_DE_TOP10_gonads_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up_gonads_1VSALL_ML, file=file.path(outputdirDE_gonads, "02_UP_prot_DE_gonads_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.up.10_gonads_1VSALL_ML, file=file.path(outputdirDE_gonads, "02_UP_prot_DE_TOP10_gonads_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


#down
write.table(DE.down_gonads_1VSALL, file=file.path(outputdirDE_gonads,"02_DOWN_prot_DE_gonads_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_gonads_1VSALL, file=file.path(outputdirDE_gonads,"02_DOWN_prot_DE_TOP10_gonads_1VSALL.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gonads_1VSALL_MG, file=file.path(outputdirDE_gonads, "02_DOWN_prot_DE_gonads_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_gonads_1VSALL_MG, file=file.path(outputdirDE_gonads, "02_DOWN_prot_DE_TOP10_gonads_1VSALL_MG.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down_gonads_1VSALL_ML, file=file.path(outputdirDE_gonads, "02_DOWN_prot_DE_gonads_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
write.table(DE.down.10_gonads_1VSALL_ML, file=file.path(outputdirDE_gonads, "02_DOWN_prot_DE_TOP10_gonads_1VSALL_ML.txt"), row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


#-------------------------------------------------------------------
#  05.  SAVING DATA
#-------------------------------------------------------------------
save(s2_gfb_counts, s2_gfb_countsF, pData, design, yF, yN, yN2, fit,
     
     contrast_gills_1VS1, qlf_gills_1VS1, resume_gills_1VS1,
     all_prot_DE_gills_1VS1, all_prot_DE_gills_1VS1_MG, all_prot_DE_gills_1VS1_ML, DE.up_gills_1VS1, DE.up_gills_1VS1_MG, DE.up_gills_1VS1_ML, DE.down_gills_1VS1, DE.down_gills_1VS1_MG, DE.down_gills_1VS1_ML,
     contrast_gills_1VSALL, qlf_gills_1VSALL, resume_gills_1VSALL,
     all_prot_DE_gills_1VSALL, all_prot_DE_gills_1VSALL_MG, all_prot_DE_gills_1VSALL_ML, DE.up_gills_1VSALL, DE.up.10_gills_1VSALL, DE.up_gills_1VSALL_MG, DE.up.10_gills_1VSALL_MG, DE.up_gills_1VSALL_ML, DE.up.10_gills_1VSALL_ML, DE.down_gills_1VSALL, DE.down.10_gills_1VSALL, DE.down_gills_1VSALL_MG, DE.down_gills_1VSALL_ML, DE.down.10_gills_1VSALL_ML,
     
     contrast_caeca_1VS1, qlf_caeca_1VS1, resume_caeca_1VS1,
     all_prot_DE_caeca_1VS1, all_prot_DE_caeca_1VS1_MG, all_prot_DE_caeca_1VS1_ML, DE.up_caeca_1VS1, DE.up_caeca_1VS1_MG, DE.up_caeca_1VS1_ML, DE.down_caeca_1VS1, DE.down_caeca_1VS1_MG, DE.down_caeca_1VS1_ML,
     contrast_caeca_1VSALL, qlf_caeca_1VSALL, resume_caeca_1VSALL,
     all_prot_DE_caeca_1VSALL, all_prot_DE_caeca_1VSALL_MG, all_prot_DE_caeca_1VSALL_ML, DE.up_caeca_1VSALL, DE.up.10_caeca_1VSALL, DE.up_caeca_1VSALL_MG, DE.up.10_caeca_1VSALL_MG, DE.up_caeca_1VSALL_ML, DE.up.10_caeca_1VSALL_ML, DE.down_caeca_1VSALL, DE.down.10_caeca_1VSALL, DE.down_caeca_1VSALL_MG, DE.down_caeca_1VSALL_ML, DE.down.10_caeca_1VSALL_ML,
     
     contrast_gonads_1VS1, qlf_gonads_1VS1, resume_gonads_1VS1,
     all_prot_DE_gonads_1VS1, all_prot_DE_gonads_1VS1_MG, all_prot_DE_gonads_1VS1_ML, DE.up_gonads_1VS1, DE.up_gonads_1VS1_MG, DE.up_gonads_1VS1_ML, DE.down_gonads_1VS1, DE.down_gonads_1VS1_MG, DE.down_gonads_1VS1_ML,
     contrast_gonads_1VSALL, qlf_gonads_1VSALL, resume_gonads_1VSALL,
     all_prot_DE_gonads_1VSALL, all_prot_DE_gonads_1VSALL_MG, all_prot_DE_gonads_1VSALL_ML, DE.up_gonads_1VSALL, DE.up.10_gonads_1VSALL, DE.up_gonads_1VSALL_MG, DE.up.10_gonads_1VSALL_MG, DE.up_gonads_1VSALL_ML, DE.up.10_gonads_1VSALL_ML, DE.down_gonads_1VSALL, DE.down.10_gonads_1VSALL, DE.down_gonads_1VSALL_MG, DE.down_gonads_1VSALL_ML, DE.down.10_gonads_1VSALL_ML,
     
     file=file.path(outputdirDE, "data_DE.Rdata"))

