setwd("Biomedical_AI_2022/R_Code")
print('Loading Required Packages')
library(ggplot2)
library(wateRmelon)
library(verification)
library(dplyr)
library(RPMM)
library(stringr)
source("functions.R")

#check versions for writeup
sessionInfo()

# Use minfi to read .idat files as RGSet
# extended RGset is required for some of the normalization methods
basedir <- '/GWAS_Source/Methylation'

print('Reading in Targets')
targets = read.csv("LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
remove <- c("LBC360558","LBC360760","LBC360536","LBC361272","LBC360213","LBC360262","LBC361264","LBC361030","LBC360412","LBC361076","LBC360721")
targets = targets[-which(targets$ID_raw %in% remove),]
targets_21 <- targets[which(targets$WAVE == 1) , ]
targets_21 <- targets_21[which(targets_21$cohort == 'LBC21') , ]
targets_36 <- targets[which(targets$WAVE == 4) , ]
targets_36 <- targets_36[which(targets_36$cohort == 'LBC36') , ]

# targets <- targets_21 #if normalizing LBC21 only
# targets <- targets_36 #if normalizing LBC36 only
targets <- rbind(targets_36,targets_21) #if normalizing both LBC21 and LBC36 combined

## The rest of the file can remain unchanged. 

print('Mapping Targets to idat files')
# Extended RGSet object
RGSet <- read.metharray.exp('Biomedical_AI_2022/idats/', targets, force=T, extended=T)
# # Non-extended object
# RGSet2 <- read.metharray.exp('Biomedical_AI_2022/idats/', targets, force=T, extended=F)

# Process the RGSet to obtain raw methylated/unmethylated intensities
MSet.raw <- mypreprocessRaw(RGSet)

#Annotation (annotate MSet using EPIC array annotation)
annot <- getAnnotation(MSet.raw)

print('Excludes SNP probes')
# Exclude SNP probes and cross-hybridising probes (McCartney et al 2016, Genomics Data)
snps <- read.table("mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", 
                   sep='\t', header=T)

cg_crosshyb <- read.table("mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt",
                          sep='\t', header=F)
ch_crosshyb <- read.table("mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt",
                          sep='\t', header=F)
crosshyb <- rbind(cg_crosshyb, ch_crosshyb)$V1 %>% as.character 
           

# Population is N.Irish - use EUR allele frequencies for snp filtering
snp_probes <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character

exclude_probes <- c(crosshyb, snp_probes) %>% unique

# Remove the SNP/Cross-hybridising probes from MSet
MyMSet <- MSet.raw[which(!rownames(MSet.raw) %in% exclude_probes), ]

# Apply pfilter 
# Probes with >1% samples with p>0.05 
# Samples with >1% probes with p > 0.05
# Sites with beadcount <3 in 5% of samples

# Which probes contain the "rs control" probes?
rs_ind <- grep("rs", rownames(MyMSet))

# get methylated intensities (rs probes not included)
mn <- minfi::getMeth(MyMSet[-rs_ind, ])
mn2 <- minfi::getMeth(MyMSet)

# get unmethylated intensities (rs probes not included)
un <- getUnmeth(MyMSet[-rs_ind, ])
un2 <- getUnmeth(MyMSet[, ])

# get beadcounts from RGSet
bc <- beadcount(RGSet)
bc <- bc[rownames(mn), ]

# Get detection p-values for each probe
detP <- detectionP(RGSet)
detP <- detP[rownames(mn), ]

# Apply p-filter to remove poor-performing probes/samples
MSet.pf <- pfilter(mn,un,bc=bc,pn=detP)

# 0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed 
# Samples removed:  
# 1365 sites were removed as beadcount <3 in 5 % of samples 
# 6759 sites having 1 % of samples with a detection p-value greater than 0.05 were removed

print(give_me_error)

keep_samps <- sampleNames(MSet.pf) # all samples are to be kept

# probes to keep along with rs probes for downstream step
keep_probes <- c(rownames(MSet.pf$mn), rownames(MyMSet)[grep("rs", rownames(MyMSet))])

# Subset mset to the probes and samples you want
MyMSet.pf <- MyMSet[keep_probes, keep_samps]

# get vector of probe types (I or II) for normalization below
onetwo <- as.character(annot[match(rownames(MyMSet.pf),rownames(annot)),"Type"])
designv <- gsub("II", "2", onetwo)
designv <- gsub("I", "1", designv)

# Get betas, use illumina offset to avoid NAs
M <- getMeth(MyMSet.pf)
U <- getUnmeth(MyMSet.pf)
betas <- M/(M+U+100)

# Here, we use a data-driven approach to determine 
# the optimum normalization function (see Pidsley et al. WateRmelon paper)

# roco variable required for some normalization methods
roco <- as.character(MyMSet.pf$Array)

print('Apply Normalisation')
# Apply different normalization functions to the data (may take a while)
Methyl.raw <- MyMSet.pf

# saveRDS(Methyl.raw, file"../LBC21_Norm_RDS_Files/Methyl21.Raw.rds")
# saveRDS(rownames(Methyl.raw), file="../LBC21_probes.rds")

# saveRDS(Methyl.raw, file"../LBC36_Norm_RDS_Files/Methyl36.Raw.rds"
# saveRDS(rownames(Methyl.raw), file="../LBC36_probes.rds")


# Single probeset required for within-array normalization consistency
probes = intersect(readRDS("../LBC21_probes.rds"), 
                             readRDS("../LBC36_probes.rds"))

# get vector of probe types (I or II) for normalization below
onetwo <- as.character(annot[match(probes,rownames(annot)),"Type"])
designv <- gsub("II", "2", onetwo)
designv <- gsub("I", "1", designv)


set.seed(1.123) # set seed as SWAN selects a random subset of probes
Methyl.swan <- mypreprocessSWAN(rgSet=RGSet[,colnames(MyMSet.pf)], mSet=MyMSet.pf[probes,], verbose=T)

Methyl.Noob <- mypreprocessNoob(RGSet[,colnames(MyMSet.pf)]) 
Methyl.Noob <- Methyl.Noob[probes,]

Methyl.Funnorm <- getBeta(preprocessFunnorm(rgSet=RGSet[,colnames(MyMSet.pf)])) # Unable to force rsids through here
Methyl.Funnorm <- rbind(Methyl.Funnorm, betas(MyMSet.pf)[grep("rs", rownames(betas(MyMSet.pf))),])[probes,] 

# At this point, I've removed RGSets to make space for the large objects generated here

Methyl.dasen <- dasen(MyMSet.pf[probes,],onetwo=onetwo)
Methyl.nasen <- nasen(MyMSet.pf[probes,])
Methyl.nanet <- nanet(MyMSet.pf[probes,])
Methyl.naten <- naten(MyMSet.pf[probes,])
Methyl.nanes <- nanes(MyMSet.pf[probes,])
Methyl.danes <- danes(MyMSet.pf[probes,])
Methyl.danet <- danet(MyMSet.pf[probes,])
Methyl.danen <- danen(MyMSet.pf[probes,],roco=pData(RGSet)[colnames(MyMSet.pf), "pos"])
Methyl.daten1 <- daten1(MyMSet.pf[probes,],roco=pData(RGSet)[colnames(MyMSet.pf), "pos"])
Methyl.daten2 <- daten2(MyMSet.pf[probes,],roco=pData(RGSet)[colnames(MyMSet.pf), "pos"])
set.seed(1.123)
Methyl.BMIQ <- BMIQ(betas(MyMSet.pf)[probes,], designv)

Methyl.fuks <- fuks(MyMSet.pf[probes,], annot[probes,])


p2 = rownames(MyMSet.pf)[which(rownames(MyMSet.pf) %in% probes)]
p2 = p2[-grep("rs", p2)]

Methyl.tost <- mynormalizeIlluminaMethylation(beta = betas(MyMSet.pf)[p2,],     # Unable to force rsids through here
                                              detect.pval = detectionP(RGSet)[p2,], 
                                              quantile.norm.pvalThreshold = NA, 
                                              probeAnnotations= annot[p2,] )

Methyl.tost$beta = rbind(Methyl.tost$beta, betas(MyMSet.pf)[grep("rs", rownames(betas(MyMSet.pf))),])[probes,]



saveRDS(MyMSet.pf[probes,], file="./../LBCComb_Norm_RDS_Files/MethylComb.Raw.rds")
saveRDS(Methyl.dasen , file = './../LBCComb_Norm_RDS_Files/MethylComb.dasen.rds')
saveRDS(Methyl.danes , file = './../LBCComb_Norm_RDS_Files/MethylComb.danes.rds')
saveRDS(Methyl.danet , file = './../LBCComb_Norm_RDS_Files/MethylComb.danet.rds')
saveRDS(Methyl.nanes , file = './../LBCComb_Norm_RDS_Files/MethylComb.nanes.rds')
saveRDS(Methyl.nanet , file = './../LBCComb_Norm_RDS_Files/MethylComb.nanet.rds')
saveRDS(Methyl.nasen , file = './../LBCComb_Norm_RDS_Files/MethylComb.nasen.rds')
saveRDS(Methyl.naten , file = './../LBCComb_Norm_RDS_Files/MethylComb.naten.rds')
saveRDS(Methyl.danen , file = './../LBCComb_Norm_RDS_Files/MethylComb.danen.rds')
saveRDS(Methyl.daten1 , file = './../LBCComb_Norm_RDS_Files/MethylComb.daten1.rds')
saveRDS(Methyl.daten2 , file = './../LBCComb_Norm_RDS_Files/MethylComb.daten2.rds')
saveRDS(Methyl.Noob , file = './../LBCComb_Norm_RDS_Files/MethylComb.noob.rds')
saveRDS(Methyl.swan , file = './../LBC36_Norm_RDS_Files/Methyl36.swan.rds')
saveRDS(Methyl.BMIQ , file = './../LBCComb_Norm_RDS_Files/MethylComb.BMIQ.rds')
saveRDS(Methyl.fuks, file="./../LBCComb_Norm_RDS_Files/MethylComb.fuks.rds")
saveRDS(Methyl.tost, file="./../LBCComb_Norm_RDS_Files/MethylComb.tost.rds")
saveRDS(Methyl.Funnorm, file="./../LBC21_Norm_RDS_Files/Methyl21.Funnorm.rds")


rm(list=ls())
gc()

set.seed(1.123); bmiq_21 <- BMIQ(betas(raw_21)[probes,], designv); set.seed(1.123); bmiq_36 <- BMIQ(betas(raw_36)[probes,], designv); set.seed(1.123); bmiq_comb <- BMIQ(betas(raw_comb)[probes,], designv); 


