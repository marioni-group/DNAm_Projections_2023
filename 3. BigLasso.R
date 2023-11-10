setwd("Daniel/Biomed_AI_Jun2022")


library(bigmemory)
library(bigmemoryExt)
library(optparse)
library(biglasso)
library(tidyverse)
library(foreign)
library(coxme)
library(kinship2)


# In case kinship matrix is to be fitted (2022-11-20)
ped <- read.csv("pedigree.csv", header=T, stringsAsFactors=F)

  ped$father <- as.numeric(ped$father)
  ped$mother <- as.numeric(ped$mother)

  ped$father[ped$father==0] <- NA
  ped$mother[ped$mother==0] <- NA

  table(ped$sex)
  ped$sex[which(ped$sex=="M")] <- "Male"
  ped$sex[which(ped$sex=="F")] <- "Female"
  

kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))

k <- kinship(kin) 


targets = readRDS("GS20k_Targets.rds")
agesex = read.csv("GS_dataset/agemonths.csv")
agesex = agesex[which(agesex$id %in% targets$Sample_Name),]
pc = read.table("GS_dataset/GS20K_ALL_MAF5_PCA.eigenvec")

bmi = read.csv("body.csv")
bmi = bmi[which(bmi$id %in% targets$Sample_Name), c("id", "bmi")]
pheno = agesex
pheno$bmi = bmi[match(agesex$id, bmi$id), "bmi"]
pheno = merge(pheno, pc, by.x="id", by.y="V2", all.x=TRUE)
pheno$bmi_res = resid(lm(log(bmi) ~ age_months + sex + V3 + V4 + V5 + V6 + 
                     V7 + V8 + V9 + V10 + V11 + V12, data=pheno, na.action=na.exclude)) 


bmi_kin = resid(lmekin(log(pheno$bmi) ~ 
             pheno$age_months + pheno$sex + pheno$V3 + pheno$V4 + pheno$V5 + pheno$V6 + 
                     pheno$V7 + pheno$V8 + pheno$V9 + pheno$V10 + pheno$V11 + pheno$V12 + (1|pheno$id), 
            varlist=k*2, na.action=na.exclude)) # %>% extract_coxme_table

pheno[as.numeric(names(bmi_kin)),"bmi_kin"] = bmi_kin

# x needs to be a list of big matrices, and they have to have the same rownames
cbindBM_list <- function(x, binding="right", 
                         z=NULL, type=NULL, separated=NULL,
                         backingfile=NULL, backingpath=NULL,
                         descriptorfile=NULL, binarydescriptor=FALSE,
                         shared=TRUE, erase = TRUE)
{
  
  if (is.null(type)) type <- typeof(x[[1]])
  if (is.big.matrix(x[[1]])) {
    if (is.null(separated)) separated <- is.separated(x[[1]])
  } else {
    separated <- FALSE
  }
  
  cols_list <- list()
  total_cols <- 0
  for (i in 1:length(x)) {
    cols <- cleanupcols(NULL, ncol(x[[i]]), colnames(x[[i]]))
    cols_list <- append(cols_list, list(cols))
    total_cols <- total_cols + ncol(x[[i]])
  }    
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x[[1]]), ncol=total_cols, type=type, init=NULL,
                    dimnames=dimnames(x[[1]]), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared=shared)
  }
  
  counter <- 0
  for (i in 1:length(cols_list)) {
    print(i)
    if (i == 1) {
      z[, 1:length(cols_list[[i]])] <- x[[i]][,cols_list[[i]]]
    } else {
      z[, (counter + 1):(counter + length(cols_list[[i]]))] <- x[[i]][,cols_list[[i]]]
    }
    counter <- counter + length(cols_list[[i]])
    print(counter)
    
    if (erase == TRUE) {
      cat("\nErasing chunk and liberating memory...\n\n")
      x[[i]] <- "Replacement"
      gc()
    }
  }
  return(z)
}
cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}
################################################################################
### PREP DNAm 
################################################################################
dat <- readRDS("/GS_20k/mvals.rds")

# NA's - mean-impute data
for(i in 1:ncol(dat)){
    if(length(which(is.infinite(dat[,i])))>=1){
    dat[which(is.infinte(dat[,i])),i] = NA
  }

  if(length(which(is.na(dat[,i])))>=1){
    dat[which(is.na(dat[,i])),i] = mean(dat[,i], na.omit=T)
  }
}

var_probes = apply(dat, 1, sd)
probes_200k = names(var_probes)[rev(order(var_probes))]

targets = readRDS("GS20k/GS20k_Targets.rds")
x = readRDS("Daniel/Biomed_AI_Jun2022/cpg_sites_Comb.rds")
y = readRDS("Daniel/Biomed_AI_Jun2022/cpg_sites_21.rds")
z = readRDS("Daniel/Biomed_AI_Jun2022/cpg_sites_36.rds")
cpg = intersect(x, intersect(y,z))

pheno$Sample_Sentrix_ID = targets[match(pheno$id, targets$Sample_Name), "Sample_Sentrix_ID"]
dat = dat[which(rownames(dat) %in% cpg),pheno$Sample_Sentrix_ID]

# 200k variable probes
probes_200k_450 = probes_200k[which(probes_200k %in% rownames(dat))][1:200000]
dat = dat[probes_200k_450, ]


bmi <- pheno
bmi = bmi[-which(is.na(bmi$bmi_res)), ]
meth <- dat[,bmi$Sample_Sentrix_ID]
div <- 15 # Number of chunks to divide OG methylation dataframe
por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()
for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}
# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)
cat("\nRAM clean up...\n\n")
gc()
cat("\nFusing chunks!\n\n")
dat2 <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)
# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(dat2)<- names
dat2 = t(dat2)
colnames(dat2)= rownames(dat)
################################################################################
### RUN BMI TRAINING
################################################################################


location <- "Daniel/Biomed_AI_Jun2022/"
set.seed(1783) # set seed to ensure fold variation minimised 
# Set list of traits to run through 
list <- "bmi_residual"

bmi$Sample_Sentrix_ID = targets[match(bmi$id, targets$Sample_Name), "Sample_Sentrix_ID"]
bmi = bmi[match(colnames(dat2), bmi$Sample_Sentrix_ID),]
identical(colnames(dat2), bmi$Sample_Sentrix_ID) # TRUE
# Assign ytrain
ytrain <- bmi
  i <- "bmi_res"  
  q <- ytrain[,"Sample_Sentrix_ID"] # Get just basenames for people in the y variable 
  p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
  name_p <- colnames(p) # Get the name of the trait for this iteration
  y <- cbind(q,p) # Bind Basename and protein data together into one set 
  names(y)[2] <- "pheno" # Assign a generic name to the protein variable
  y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model
  # Run training 
  lasso.cv <- cv.biglasso(dat2, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
  fit <- biglasso(dat2, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 
  coefs <- coef(fit) # Extract coeficients 
  coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  coefs$Predictor <- name_p # Assign protein identifier
  names(coefs)[1] <- "Coefficient" # Tidy naming 
  coefs$CpG <- rownames(coefs) # Create episcores column
  # coefs <- coefs[-1,]
  write.csv(coefs, file = paste0(location, name_p, "GS20k_train_weights_bmi_200k.csv"), row.names = F)

# Assign ytrain
ytrain <- bmi
  i <- "bmi_kin"  
  q <- ytrain[,"Sample_Sentrix_ID"] # Get just basenames for people in the y variable 
  p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
  name_p <- colnames(p) # Get the name of the trait for this iteration
  y <- cbind(q,p) # Bind Basename and protein data together into one set 
  names(y)[2] <- "pheno" # Assign a generic name to the protein variable
  y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model
  # Run training 
  lasso.cv <- cv.biglasso(dat2, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
  fit <- biglasso(dat2, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 
  coefs <- coef(fit) # Extract coeficients 
  coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  coefs$Predictor <- name_p # Assign protein identifier
  names(coefs)[1] <- "Coefficient" # Tidy naming 
  coefs$CpG <- rownames(coefs) # Create episcores column
  # coefs <- coefs[-1,]
  write.csv(coefs, file = paste0(location, name_p, "GS20k_train_weights_bmi_kinship_200k.csv"), row.names = F)

