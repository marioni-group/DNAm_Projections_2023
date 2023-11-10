setwd("Biomedical_AI_2022/R_Code")

library(ggplot2)
library(wateRmelon)
library(verification)
library(dplyr)
library(RPMM)
library(stringr)
source("functions.R")
#check versions for writeup
sessionInfo()

#Annotation (annotate MSet using EPIC array annotation)
annot <- readRDS("Daniel/450k_annotation.rds")
targets = read.csv("LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")



	# Apply different normalization functions to the data (may take a while)
	Methyl.raw <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.Raw.rds")
	Methyl.swan <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.swan.rds")
	Methyl.Noob <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.noob.rds")
	Methyl.Funnorm <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.Funnorm.rds")
	Methyl.dasen <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.dasen.rds")
	Methyl.nasen <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.nasen.rds")
	Methyl.nanet <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.nanet.rds")
	Methyl.naten <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.naten.rds")
	Methyl.nanes <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.nanes.rds")
	Methyl.danes <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.danes.rds")
	Methyl.danet <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.danet.rds")
	Methyl.danen <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.danen.rds")
	Methyl.daten1 <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.daten1.rds")
	Methyl.daten2 <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.daten2.rds")
	Methyl.BMIQ <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.BMIQ.rds")
	Methyl.fuks <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.fuks.rds")
	Methyl.tost <- readRDS("../LBCComb_Norm_RDS_Files/MethylComb.tost.rds")


#DMRSE: DMRSE Metric (see Pidsley et al)
DMRSE_raw <- dmrse_row(Methyl.raw)
DMRSE_swan <- dmrse_row(Methyl.swan)
DMRSE_Noob <- dmrse_row(Methyl.Noob)
DMRSE_Funnorm <- dmrse_row(Methyl.Funnorm)
DMRSE_dasen <- dmrse_row(Methyl.dasen)
DMRSE_nasen <- dmrse_row(Methyl.nasen)
DMRSE_nanet <- dmrse_row(Methyl.nanet)
DMRSE_naten <- dmrse_row(Methyl.naten)
DMRSE_nanes <- dmrse_row(Methyl.nanes)
DMRSE_danes <- dmrse_row(Methyl.danes)
DMRSE_danet <- dmrse_row(Methyl.danet)
DMRSE_danen <- dmrse_row(Methyl.danen)
DMRSE_daten1 <- dmrse_row(Methyl.daten1)
DMRSE_daten2 <- dmrse_row(Methyl.daten2)
DMRSE_fuks <- dmrse_row(Methyl.fuks)
DMRSE_BMIQ <- dmrse_row(Methyl.BMIQ)
DMRSE_tost <- dmrse_row(Methyl.tost$beta)

#GCOSE: GCOSE Metric (see Pidsley et al)
GCOSE_raw <- genki(Methyl.raw,se=T)
GCOSE_swan <- genki(Methyl.swan, se=T)
GCOSE_Noob <- genki(Methyl.Noob, se=T)
GCOSE_Funnorm <- genki(Methyl.Funnorm, se=T)
GCOSE_dasen <- genki(Methyl.dasen, se=T)
GCOSE_nasen <- genki(Methyl.nasen, se=T)
GCOSE_nanet <- genki(Methyl.nanet, se=T)
GCOSE_naten <- genki(Methyl.naten, se=T)
GCOSE_nanes <- genki(Methyl.nanes, se=T)
GCOSE_danes <- genki(Methyl.danes, se=T)
GCOSE_danet <- genki(Methyl.danet, se=T)
GCOSE_danen <- genki(Methyl.danen, se=T)
GCOSE_daten1 <- genki(Methyl.daten1, se=T)
GCOSE_daten2 <- genki(Methyl.daten2, se=T)
GCOSE_fuks <- genki(Methyl.fuks, se=T)
GCOSE_BMIQ <- genki(Methyl.BMIQ, se=T)

# GCOSE_tost <- genki(rbind(Methyl.tost$beta, betas(Methyl.raw)[grep("rs", rownames(Methyl.raw)),]), se=T)# Append raw rs probes to tost (lbc21/36)
GCOSE_tost <- genki(Methyl.tost$beta, se=T)# Append raw rs probes to tost (combined)


#Seabird metric (see Pidsley et al)
bnraw <- betas(Methyl.raw)
bnswan <- betas(Methyl.swan)
bnNoob <- betas(Methyl.Noob)
bnFunnorm <- Methyl.Funnorm
bndasen <- betas(Methyl.dasen)
bnnasen <- betas(Methyl.nasen)
bnnanet <- betas(Methyl.nanet)
bnnaten <- betas(Methyl.naten)
bnnanes <- betas(Methyl.nanes)
bndanes <- betas(Methyl.danes)
bndanet <- betas(Methyl.danet)
bndanen <- betas(Methyl.danen)
bndaten1 <- betas(Methyl.daten1)
bndaten2 <- betas(Methyl.daten2)
bnfuks <- Methyl.fuks
bnBMIQ <- Methyl.BMIQ
bntost <- Methyl.tost$beta

# xchr <- rownames(Methyl.raw) %in% rownames(annot[which(annot$chr=="chrX"),])
sex = targets[match(colnames(Methyl.raw), targets$Basename),"sex"]
sex = gsub("F", 2, sex)
sex = gsub("M", 1, sex)
sex = as.numeric(sex) 
# one missing sex: imputed as male based on mean Xchr value
sex[is.na(sex)] = 1

SEABIRD_raw <- mySeabi(bnraw, sex=sex, X=rownames(bnraw) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_swan<- mySeabi(bnswan, sex=sex, X=rownames(bnswan) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_Noob<- mySeabi(bnNoob, sex=sex, X=rownames(bnNoob) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_Funnorm<- mySeabi(bnFunnorm, sex=sex, X=rownames(bnFunnorm) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_dasen <- mySeabi(bndasen, sex=sex, X=rownames(bndasen) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_nasen <- mySeabi(bnnasen, sex=sex, X=rownames(bnnasen) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_nanet <- mySeabi(bnnanet, sex=sex, X=rownames(bnnanet) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_naten <- mySeabi(bnnaten, sex=sex, X=rownames(bnnaten) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_nanes <- mySeabi(bnnanes, sex=sex, X=rownames(bnnanes) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_danes <- mySeabi(bndanes, sex=sex, X=rownames(bndanes) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_danet <- mySeabi(bndanet, sex=sex, X=rownames(bndanet) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_danen <- mySeabi(bndanen, sex=sex, X=rownames(bndanen) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_daten1 <- mySeabi(bndaten1, sex=sex, X=rownames(bndaten1) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_daten2 <- mySeabi(bndaten2, sex=sex, X=rownames(bndaten2) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_fuks <- mySeabi(bnfuks, sex=sex, X=rownames(bnfuks) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_BMIQ <- mySeabi(bnBMIQ, sex=sex, X=rownames(bnBMIQ) %in% rownames(annot[which(annot$chr=="chrX"),]))
SEABIRD_tost <- mySeabi(bntost, sex=sex, X=rownames(bntost) %in% rownames(annot[which(annot$chr=="chrX"),]))

# Tabulate performance of each metric to determine best normalization method for your data

DMRSE <- c(DMRSE_raw, DMRSE_swan, DMRSE_Noob, DMRSE_Funnorm, DMRSE_dasen,
         DMRSE_nasen, DMRSE_nanet, DMRSE_naten, DMRSE_nanes, DMRSE_danes,
         DMRSE_danet, DMRSE_danen, DMRSE_daten1, DMRSE_daten2, DMRSE_fuks, DMRSE_BMIQ, DMRSE_tost)

GCOSE <- c(mean(GCOSE_raw), mean(GCOSE_swan), mean(GCOSE_Noob), mean(GCOSE_Funnorm),
         mean(GCOSE_dasen), mean(GCOSE_nasen), 
		 mean(GCOSE_nanet), mean(GCOSE_naten), mean(GCOSE_nanes), 
		 mean(GCOSE_danes), mean(GCOSE_danet), mean(GCOSE_danen), 
		 mean(GCOSE_daten1), mean(GCOSE_daten2), mean(GCOSE_fuks), mean(GCOSE_BMIQ), mean(GCOSE_tost))

SEABIRD<-c(SEABIRD_raw, SEABIRD_swan, SEABIRD_Noob, SEABIRD_Funnorm, SEABIRD_dasen, 
           SEABIRD_nasen, SEABIRD_nanet, SEABIRD_naten, SEABIRD_nanes, 
		   SEABIRD_danes, SEABIRD_danet, SEABIRD_danen, SEABIRD_daten1, 
		   SEABIRD_daten2, SEABIRD_fuks, SEABIRD_BMIQ, SEABIRD_tost)

Normalization<-c('Raw','SWAN','Noob', 'Funnorm', 'dasen','nasen','nanet','naten','nanes','danes','danet','danen','daten1','daten2', 'fuks', 'BMIQ', 'tost')
Table1<-data.frame(Normalization,"Rank of DMRSE"=rank(DMRSE),
"Rank of GCOSE"=rank(GCOSE), "Rank of Seabird"=rank(SEABIRD))
Table2<-data.frame(Normalization,Table1[,2:4],"Mean of Ranks"=rowMeans(
Table1[,2:4]), "Ranked Means"=rank(rowMeans(Table1[,2:4])))

# Save table of metric rankings. 
# "Ranked Means" column is used to identify the best normalization method

  # write.table(Table2, file="Normalisation_performance_LBC1921_PostReview.txt", sep='\t', row.names=F)
# write.table(Table2, file="Normalisation_performance_LBC1936_PostReview.txt", sep='\t', row.names=F)
 write.table(Table2, file="Normalisation_performance_LBCCombined_PostReview.txt", sep='\t', row.names=F)

# For Supplementary plot
pd_performance = data.frame(Method = Normalization, DMRSE = DMRSE, GCOSE = GCOSE, SEABIRD = SEABIRD)
# write.table(pd_performance, file="Normalisation_performance_plotdata_LBC1921_PostReview.txt", sep='\t', row.names=F)
# write.table(pd_performance, file="Normalisation_performance_plotdata_LBC1936_PostReview.txt", sep='\t', row.names=F)
 write.table(pd_performance, file="Normalisation_performance_plotdata_LBCCombined_PostReview.txt", sep='\t', row.names=F)