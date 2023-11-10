.isRG <- function(object) {
    if(!is(object, "RGChannelSet"))
        stop(sprintf("object is of class '%s', but needs to be of class 'RGChannelSet' or 'RGChannelSetExtended'",
                     class(object)))
}

.isRGOrStop <- function(object) {
    if(!is(object, "RGChannelSet"))
        stop(sprintf("object is of class '%s', but needs to be of class 'RGChannelSet' or 'RGChannelSetExtended'",
                     class(object)))
}

.getManifestString <- function(annotation) {
    if(length(annotation) == 1) {
        if(annotation == "Unknown")
            stop("Cannot get Manifest object for an 'Unknown' array")
        return(paste0(annotation, "manifest"))
    }
    if("array" %in% names(annotation)) {
        if(annotation["array"] == "Unknown")
            stop("Cannot get Manifest object for an 'Unknown' array")
        return(paste0(annotation["array"], "manifest"))
    }
    stop("unable to get the manifest string for this object")
}

  mypreprocessRaw <- 
function (rgSet) 
{
    .isRGOrStop(rgSet)
    
    probeinfo_i <- getProbeInfo(rgSet, type = "I")
    probeinfo_ii <- getProbeInfo(rgSet, type = "II")
    probeinfo_snpi <- getProbeInfo(rgSet, type = "SnpI")
    probeinfo_snpii <- getProbeInfo(rgSet, type = "SnpII")
    probeinfo_igreen <- getProbeInfo(rgSet, type = "I-Green")
    probeinfo_ired <- getProbeInfo(rgSet, type = "I-Red")

    locusNames <- c(probeinfo_i$Name, probeinfo_ii$Name, 
                       probeinfo_snpi$Name, probeinfo_snpii$Name)
                       
    M <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames), 
        dimnames = list(locusNames, colnames(rgSet)))
    U <- M
    TypeII.Name <- getProbeInfo(rgSet, type = "II")$Name
    TypeIISNP.Name <- getProbeInfo(rgSet, type = "SnpII")$Name
    M[TypeII.Name, ] <- getGreen(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA, 
        ]
    M[TypeIISNP.Name, ] <- getGreen(rgSet)[getProbeInfo(rgSet, type = "SnpII")$AddressA,
        ]
    U[TypeII.Name, ] <- getRed(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA, 
        ]
    U[TypeIISNP.Name, ] <- getRed(rgSet)[getProbeInfo(rgSet, type = "SnpII")$AddressA,
        ]
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeI.SNP <- getProbeInfo(rgSet, type = "SnpI")
    M[TypeI.Red$Name, ] <- getRed(rgSet)[TypeI.Red$AddressB, 
        ]
    M[TypeI.Green$Name, ] <- getGreen(rgSet)[TypeI.Green$AddressB, 
        ]
    M[TypeI.SNP$Name,] <- getGreen(rgSet)[TypeI.SNP$AddressB,
        ]
    U[TypeI.Red$Name, ] <- getRed(rgSet)[TypeI.Red$AddressA, 
        ]
    U[TypeI.Green$Name, ] <- getGreen(rgSet)[TypeI.Green$AddressA, 
        ]
    U[TypeI.SNP$Name, ] <- getGreen(rgSet)[TypeI.SNP$AddressA,
        ]
    out <- MethylSet(Meth = M, Unmeth = U, colData = colData(rgSet), 
        annotation = annotation(rgSet), metadata = metadata(rgSet))
    out@preprocessMethod <- c(rg.norm = "Raw (no normalization or bg correction)", 
        minfi = as.character(packageVersion("minfi")), manifest = as.character(packageVersion(.getManifestString(rgSet@annotation))))
    out
}
 



bgIntensitySwan <- function (rgSet)
{
    grnMed <- colMedians(getGreen(rgSet)[getControlAddress(rgSet,
        controlType = "NEGATIVE"), ])
    redMed <- colMedians(getRed(rgSet)[getControlAddress(rgSet,
        controlType = "NEGATIVE"), ])
    return(rowMeans(cbind(grnMed, redMed)))
}

getSubset <- function(counts, subset){
    x <- numeric(0)
    for(i in 1:3){
        x <- c(x,sample(seq(1, length(counts), by = 1)[counts == i], subset))
    }
    return(seq(1, length(counts)) %in% x)
}

normaliseChannel <- function(intensityI, intensityII, xNormSet, bg) {
    xTarget <- aveQuantile(list(intensityI[xNormSet[[1]]], intensityII[xNormSet[[2]]]))
    xNorm <- unlist(subsetQuantileNorm(list(intensityI, intensityII), xNormSet, xTarget, bg))
    names(xNorm) <- c(names(intensityI), names(intensityII))
    xNorm
}

aveQuantile <- function(X) {
    nbrOfChannels <- length(X)
    if (nbrOfChannels == 1) {
        return(X)
    }
    nbrOfObservations <- unlist(lapply(X, FUN = length), use.names = FALSE)
    maxNbrOfObservations <- max(nbrOfObservations)
    if (maxNbrOfObservations == 1) {
        return(X)
    }
    ## nbrOfFiniteObservations <- rep(maxNbrOfObservations, times = nbrOfChannels)
    quantiles <- (0:(maxNbrOfObservations - 1))/(maxNbrOfObservations - 1)
    xTarget <- vector("double", maxNbrOfObservations)
    for (cc in 1:nbrOfChannels) {
        Xcc <- X[[cc]]
        Scc <- sort(Xcc)
        nobs <- length(Scc)
        if (nobs < maxNbrOfObservations) {
            ## tt <- !is.na(Xcc)
            bins <- (0:(nobs - 1))/(nobs - 1)
            Scc <- approx(x = bins, y = Scc, xout = quantiles,ties = "ordered")$y
        }
        xTarget <- xTarget + Scc
    }
    rm(Scc, Xcc)
    xTarget <- xTarget/nbrOfChannels
    xTarget
}

subsetQuantileNorm <- function(x, xNormSet, xTarget, bg) {
    for(i in 1:length(x)){
        n <- length(x[[i]])
        nTarget <- length(xTarget)
        nNormSet <- sum(xNormSet[[i]])
        
        if(nNormSet != nTarget){
            targetQuantiles <- (0:(nTarget - 1))/(nTarget - 1)
            r <- rank(x[xNormSet[,i], i])
            xNew <-(r - 1)/(nNormSet - 1)
            xNew <- xNew[order(xNew)]
            xNorm <- approx(x = targetQuantiles, y = xTarget, xout = xNew, ties = "ordered", rule = 2)$y
        } else {
            xNorm<-xTarget
        }
        
        r <- rank(x[[i]])
        xNew <-(r - 1)/(n - 1)
        quantiles <- xNew[xNormSet[[i]]]
        quantiles <- quantiles[order(quantiles)]
        xmin <- min(x[[i]][xNormSet[[i]]]) #get min value from subset
        xmax <- max(x[[i]][xNormSet[[i]]]) #get max value from subset
        kmax <- which(xNew > max(quantiles)) 
        kmin<- which(xNew < min(quantiles))
        offsets.max <- x[[i]][kmax]-xmax
        offsets.min <- x[[i]][kmin]-xmin
        x[[i]] <- approx(x = quantiles, y = xNorm, xout = xNew, ties = "ordered")$y #interpolate
        x[[i]][kmax]<- max(xNorm) + offsets.max
        x[[i]][kmin]<- min(xNorm) + offsets.min
        x[[i]] = ifelse(x[[i]] <= 0, bg, x[[i]])    
    }  
    x
}

 mypreprocessSWAN <-
function (rgSet, mSet = NULL, verbose = FALSE) 
{
    if (is.null(mSet)) 
        mSet <- mypreprocessRaw(rgSet)
    typeI <- getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")]
    typeISnp <- getProbeInfo(rgSet, type = "SnpI")[, c("Name", "nCpG")]
    typeII <- getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")]
    typeIISnp <- getProbeInfo(rgSet, type = "SnpII")[, c("Name", "nCpG")]
    CpG.counts <- rbind(typeI, typeII, typeISnp, typeIISnp)
    CpG.counts$Name <- as.character(CpG.counts$Name)
    CpG.counts$Type <- rep(c("I", "I", "II", "II"), times = c(nrow(typeI), nrow(typeISnp),
        nrow(typeII), nrow(typeIISnp)))
    names(CpG.counts)[2] <- "CpGs"
    counts <- CpG.counts[CpG.counts$Name %in% rownames(mSet), 
        ]
    subset <- min(table(counts$CpGs[counts$Type == "I" & counts$CpGs %in% 
        1:3]), table(counts$CpGs[counts$Type == "II" & counts$CpGs %in% 
        1:3]))
    bg <- bgIntensitySwan(rgSet)
    methData <- getMeth(mSet)
    unmethData <- getUnmeth(mSet)
    xNormSet <- vector("list", 2)
    xNormSet[[1]] <- getSubset(counts$CpGs[counts$Type == "I"], 
        subset)
    xNormSet[[2]] <- getSubset(counts$CpGs[counts$Type == "II"], 
        subset)
    normMethData <- matrix(NA_real_, ncol = ncol(methData), nrow = nrow(methData))
    colnames(normMethData) <- colnames(mSet)
    normUnmethData <- normMethData
    normSet <- mSet
    for (i in 1:ncol(mSet)) {
        if (verbose) 
            message(sprintf("[preprocessSwan] Normalizing array %d of %d\n", 
                i, ncol(mSet)))
        normMeth <- normaliseChannel(methData[rownames(methData) %in% 
            counts$Name[counts$Type == "I"], i], methData[rownames(methData) %in% 
            counts$Name[counts$Type == "II"], i], xNormSet, bg[i])
        normMethData[, i] <- normMeth
        normUnmeth <- normaliseChannel(unmethData[rownames(unmethData) %in% 
            counts$Name[counts$Type == "I"], i], unmethData[rownames(unmethData) %in% 
            counts$Name[counts$Type == "II"], i], xNormSet, bg[i])
        normUnmethData[, i] <- normUnmeth
    }
    rownames(normMethData) <- names(normMeth)
    rownames(normUnmethData) <- names(normUnmeth)
    assay(normSet, "Meth") <- normMethData
    assay(normSet, "Unmeth") <- normUnmethData
    normSet@preprocessMethod <- c(rg.norm = sprintf("SWAN (based on a MethylSet preprocesses as '%s'", 
        preprocessMethod(mSet)[1]), minfi = as.character(packageVersion("minfi")), 
        manifest = as.character(packageVersion(.getManifestString(annotation(rgSet)))))
    normSet
}



mypreprocessNoob <- function(rgSet, offset=15, dyeCorr=TRUE, verbose = TRUE,
                           dyeMethod = c("single", "reference")) {
    .isRGOrStop(rgSet)
    subverbose <- max(as.integer(verbose) - 1L, 0)
    dyeMethod <- match.arg(dyeMethod)

    ## Extraction of the out-of-band controls
    controls <- getOOB(rgSet)
    names(controls) <- c("Cy3", "Cy5")
    mset <- mypreprocessRaw(rgSet)
    meth <- getMeth(mset)
    unmeth <- getUnmeth(mset)

    if (any(meth<=0)){
        meth[which(meth<=0)] <- 1
    }
    if (any(unmeth<=0)){
        unmeth[which(unmeth<=0)] <- 1
    }

    probe.type <- getProbeType(mset, withColor=TRUE)
    cy3.probes <- which(probe.type=="IGrn")
    cy5.probes <- which(probe.type=="IRed")
    d2.probes <- which(probe.type=="II")

    dat <- list(Cy3 = list(M =  as.matrix(meth[cy3.probes,]), 
                           U =  as.matrix(unmeth[cy3.probes,]),
                           D2 = as.matrix(meth[d2.probes,])), 
                Cy5 = list(M =  as.matrix(meth[cy5.probes,]), 
                           U =  as.matrix(unmeth[cy5.probes,]),
                           D2 = as.matrix(unmeth[d2.probes,])))

    rows <- lapply(dat, function(ch) {
        sapply(names(ch), function(nch) {
            nrow(ch[[nch]])
        })
    })
    last <- lapply(rows, cumsum)
    first <- lapply(names(last), function(nch) {
        last[[nch]] - rows[[nch]] + 1
    })
    names(first) <- names(last)

    estimates <- lapply(names(dat), function(nch) { 
        xf <- rbind(dat[[nch]][['M']], dat[[nch]][['U']], dat[[nch]][['D2']])
        xs <- normexp.get.xs(xf = xf, controls = controls[[nch]],
                             offset=offset, verbose = subverbose)
        names(xs[['params']]) <- paste(names(xs[['params']]), nch, sep='.')
        names(xs[['meta']]) <- paste(names(xs[['meta']]), nch, sep='.')
        xs
    })
    names(estimates) <- names(dat)

    if (length(cy3.probes)>0){
        cy3.M <- first[['Cy3']][['M']]:last[['Cy3']][['M']]
        meth[cy3.probes, ] <- estimates[['Cy3']][['xs']][cy3.M,]
        cy3.U <- first[['Cy3']][['U']]:last[['Cy3']][['U']]
        unmeth[cy3.probes,] <- estimates[['Cy3']][['xs']][cy3.U,]
    }

    if (length(cy5.probes)>0){
        cy5.M <- first[['Cy5']][['M']]:last[['Cy5']][['M']]
        meth[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.M,]
        cy5.U <- first[['Cy5']][['U']]:last[['Cy5']][['U']]
        unmeth[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.U,]
    }

    if (length(d2.probes)>0){
        d2.M <- first[['Cy3']][['D2']]:last[['Cy3']][['D2']]
        d2.U <- first[['Cy5']][['D2']]:last[['Cy5']][['D2']]
        meth[d2.probes,] <- estimates[['Cy3']][['xs']][d2.M,]
        unmeth[d2.probes,] <- estimates[['Cy5']][['xs']][d2.U,]
    }
    ## This next code block does nothing because the rgSet is not returned
    ## and colData(rgSet) is not referenced below
    ## FIXME: either remove or modify the return value
    for(ch in names(estimates)) { 
        chnames <- names(estimates[[ch]][['params']])
        for(nm in chnames)
            colData(rgSet)[,nm] <- estimates[[ch]][['params']][[nm]]
    } 

    ## Performing dye bias normalization
    ## 
    ## "single" = just reciprocate out the dye bias, don't use a reference.
    ##            (similar to, but implemented differently from, unmaintained
    ##             "asmn" package by Decker et al., doi:10.4161/epi.26037)
    ## "reference" = use the least-worst sample in the batch (previous default) 
    ## 
    ## "single" is now the default: it provides single-sample preprocessing
    ## and betas/M-values produced by this method are identical to those from
    ## the "reference" version used in (e.g.) the TCGA data processing pipeline.
    ## 
    ## --tjt, 2016-06-16
    ## 
    if (dyeCorr){
        ## Background correct the Illumina normalization controls:
        ctrls <- getProbeInfo(rgSet, type = "Control")
        ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
        redControls <- getRed(rgSet)[ctrls$Address,,drop=FALSE]
        greenControls <- getGreen(rgSet)[ctrls$Address,,drop=FALSE]
        rownames(redControls) <- rownames(greenControls) <- ctrls$Type
        internal.controls <- list(Cy3 = greenControls, Cy5 = redControls)
        xcs <- lapply(names(internal.controls), function(nch) {
          xcf <- as.matrix(internal.controls[[nch]])
          normexp.get.xcs(xcf = xcf, params=estimates[[nch]][['params']])
        })
        names(xcs) <- names(dat)
        internal.controls[['Cy3']] <- xcs[["Cy3"]]
        internal.controls[['Cy5']] <- xcs[["Cy5"]]

        if (rgSet@annotation[["array"]]=="IlluminaHumanMethylation450k" || 
                rgSet@annotation[["array"]]=="IlluminaHumanMethylationEPIC"){
            CG.controls <- rownames(internal.controls[[1]]) %in% c("NORM_C", "NORM_G")
            AT.controls <- rownames(internal.controls[[1]]) %in% c("NORM_A", "NORM_T")
        } else {
            CG.controls <- rownames(internal.controls[[1]]) %in% c("Normalization-Green")
            AT.controls <- rownames(internal.controls[[1]]) %in% c("Normalization-Red")
        }

        ## Dye bias normalization with the corrected Illumina control probes:
        Green.avg <- colMeans(internal.controls[["Cy3"]][CG.controls,, drop=FALSE])
        Red.avg <- colMeans(internal.controls[["Cy5"]][AT.controls,, drop=FALSE])
        R.G.ratio <- Red.avg/Green.avg

        if (dyeMethod == "single") {
          if(verbose) {
            cat('[MypreprocessNoob] Applying R/G ratio flip to fix dye bias...\n')
          }
          Red.factor <- 1 / R.G.ratio
          Grn.factor <- 1
        } else if(dyeMethod == "reference") {
          reference <- which.min(abs(R.G.ratio-1) )
          if(verbose) {
            cat('[mypreprocessNoob] Using sample number', reference, 
                'as reference level...\n')
          }
          ref <- (Green.avg + Red.avg)[reference]/2
          if(is.na(ref)) {
              stop("'reference' refers to an array that is not present")
          }
          Grn.factor <- ref/Green.avg
          Red.factor <- ref/Red.avg
        } else { stop("unknown 'dyeMethod'") }

        Grn <- list(M = as.matrix(meth[cy3.probes,]), 
                    U = as.matrix(unmeth[cy3.probes,]),
                    D2 = as.matrix(meth[d2.probes,]))
        Red <- list(M = as.matrix(meth[cy5.probes,]), 
                    U = as.matrix(unmeth[cy5.probes,]),
                    D2 = as.matrix(unmeth[d2.probes,]))

        ## do this regardless of reference or equalization approach
        Red <- lapply(Red, function(y) sweep(y, 2, FUN="*", Red.factor))
        meth[cy5.probes,] <- Red$M
        unmeth[cy5.probes,] <- Red$U
        unmeth[d2.probes,] <- Red$D2

        ## but only adjust the green channel if using the older reference method
        if (dyeMethod == "reference") {
          Grn <- lapply(Grn, function(y) sweep(y, 2, FUN="*", Grn.factor))
          meth[cy3.probes,] <- Grn$M
          unmeth[cy3.probes,] <- Grn$U
          meth[d2.probes,] <- Grn$D2
        }
    }

    assay(mset, "Meth") <- meth
    assay(mset, "Unmeth") <- unmeth

    mset@preprocessMethod <- c( mu.norm = 
                                    sprintf("Noob, dyeCorr=%s, dyeMethod=%s", 
                                            dyeCorr, dyeMethod))
    return(mset)
}

detachAllPackages <- function() {

  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")

  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]

  package.list <- setdiff(package.list,basic.packages)

  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}


normexp.get.xs <- function(xf, controls, offset=50, verbose = FALSE){
    if(verbose)
        message("[normexp.get.xs] Background mean & SD estimated from", nrow(controls), "probes\n")
    mu <- sigma <- alpha <- rep(NA, ncol(xf))
    for( i in 1:ncol(xf) ) {
        ests <- huber(controls[, i]) # from MASS
        mu[i] <- ests$mu
        sigma[i] <- ests$s
        alpha[i] <- max(huber(xf[, i])$mu - mu[i], 10)
    }
    pars <- data.frame(mu=mu, lsigma=log(sigma), lalpha=log(alpha))
    for(i in 1:ncol(xf))
        xf[,i] <- normexp.signal(as.numeric(pars[i,]), xf[,i]) # from limma
    return(list(xs=xf+offset, 
                params=data.frame(mu=mu, sigma=sigma, alpha=alpha, offset=offset),
                meta=c('background mean','background SD','signal mean','offset')))
}

normexp.get.xcs <- function(xcf, params){
    stopifnot(any(grepl("mu", names(params))))
    stopifnot(any(grepl("sigma", names(params))))
    stopifnot(any(grepl("alpha", names(params))))
    stopifnot(any(grepl("offset", names(params))))
    pars <- data.frame(mu=params[[grep("mu", names(params), value=TRUE)]],
                       sigma = log(params[[grep("sigma", names(params), value=TRUE)]]),
                       alpha = log(params[[grep("alpha", names(params), value=TRUE)]]))
    for(i in 1:ncol(xcf))
        xcf[,i] <- normexp.signal(as.numeric((pars[i,])), xcf[,i] ) # from limma
    return( xcf + params[[grep('offset', names(params), value=TRUE)]][1] )
}





DoPBC <-
function(beta.m,design.v){

  mval.m <- log2(beta.m/(1-beta.m));
  type1.idx <- which(design.v==1)
  type2.idx <- which(design.v==2)  
  mvalT.m <- mval.m;
  for(s in 1:ncol(beta.m)){

    neg.idx <- which(mval.m[,s]<0);
    pos.idx <- which(mval.m[,s]>0);

    neg1.idx <- intersect(neg.idx,type1.idx);
    neg2.idx <- intersect(neg.idx,type2.idx);

    pos1.idx <- intersect(pos.idx,type1.idx);
    pos2.idx <- intersect(pos.idx,type2.idx);

    d.o <- density(mval.m[neg2.idx,s],kernel="gaussian",bw=0.5);
    peakU2 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[pos2.idx,s],kernel="gaussian",bw=0.5);
    peakM2 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[neg1.idx,s],kernel="gaussian",bw=0.5);
    peakU1 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[pos1.idx,s],kernel="gaussian",bw=0.5);
    peakM1 <- abs(d.o$x[which.max(d.o$y)]);

    mvalT.m[neg2.idx,s] <- (mval.m[neg2.idx,s]/peakU2)*peakU1;
    mvalT.m[pos2.idx,s] <- (mval.m[pos2.idx,s]/peakM2)*peakM1;
    print(paste("Done for sample ",s,sep=""));
  }

  betaT.m <- 2^mvalT.m/(2^mvalT.m+1);
  return(betaT.m);
}

rowTTestPVal <- function(x,y) {
  # Based on the Welch's two test implemented by t.test.default
  # Extended to take x and y as matrices with different datasets in rows
  if (dim(x)[1] != dim(y)[1]) {
    stop("x and y must have equal numbers of rows.")
  }
  nx <- dim(x)[2]
  mx <- rowMeans(x)
  vx <- (rowSums(x^2)-mx^2*nx)/(nx-1)
  ny <- dim(y)[2]
  my <- rowMeans(y)
  vy <- (rowSums(y^2)-my^2*ny)/(ny-1)
  stderrx <- sqrt(vx/nx)
  stderry <- sqrt(vy/ny)
  stderr <- sqrt(stderrx^2 + stderry^2)
  df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 
                                                    1))
  tstat <- (mx - my)/stderr
  pval <- 2 * pt(-abs(tstat), df)
}

mySextest <- function (betas, sex,  ...) 
{
  sex <- as.integer(factor(sex))
  rowTTestPVal(betas[,sex==1],betas[,sex==2])
}

mySeabi <- function(bn, stop = 1, sex, X) {
  pVals <- mySextest(bn, sex)
  if (stop != 1) {
    1 - seabird(pVals, stop, X)
  } else {
    verification::roc.area(X,pVals)$A
  }
}
setGeneric("mySeabi")
setMethod("mySeabi","MethyLumiSet",function(bn, stop = 1, sex, X) {
  mySeabi(betas(bn),stop,sex,X)
})

mypreprocessFunnorm <- 
function (rgSet, nPCs = 2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, 
    keepCN = TRUE, ratioConvert = TRUE, verbose = TRUE) 
{
    .isRGOrStop(rgSet)
    rgSet <- updateObject(rgSet)
    if (bgCorr) {
        if (verbose && dyeCorr) {
            message("[preprocessFunnorm] Background and dye bias correction with noob")
        }
        else {
            message("[preprocessFunnorm] Background correction with noob")
        }
        gmSet <- mypreprocessNoob(rgSet, dyeCorr = dyeCorr)
        if (verbose) 
            message("[preprocessFunnorm] Mapping to genome")
        gmSet <- mapToGenome(gmSet)
    }
    else {
        if (verbose) 
            message("[preprocessFunnorm] Mapping to genome")
        gmSet <- mapToGenome(rgSet)
    }
    subverbose <- max(as.integer(verbose) - 1L, 0)
    if (verbose) 
        message("[preprocessFunnorm] Quantile extraction")
    extractedData <- .extractFromRGSet450k(rgSet)
    rm(rgSet)
    if (is.null(sex)) {
        gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
        sex <- rep(1L, length(gmSet$predictedSex))
        sex[gmSet$predictedSex == "F"] <- 2L
    }
    if (verbose) 
        message("[preprocessFunnorm] Normalization")
    if (keepCN) {
        CN <- getCN(gmSet)
    }
    gmSet <- .normalizeFunnorm450k(object = gmSet, extractedData = extractedData, 
        sex = sex, nPCs = nPCs, verbose = subverbose)
    preprocessMethod <- c(preprocessMethod(gmSet), mu.norm = sprintf("Funnorm, nPCs=%s", 
        nPCs))
    if (ratioConvert) {
        grSet <- ratioConvert(gmSet, type = "Illumina", keepCN = keepCN)
        if (keepCN) {
            assay(grSet, "CN") <- CN
        }
        grSet@preprocessMethod <- preprocessMethod
        return(grSet)
    }
    else {
        gmSet@preprocessMethod <- preprocessMethod
        return(gmSet)
    }
}


# Modified tost norm method to bypass limitations with newer minfi versions
mynormalizeIlluminaMethylation = function (beta, detect.pval, quantile.norm.pvalThreshold = 0.01,
    probeAnnotations, annotCategory = "relationToCpG")
{
    infiniumI <- probeAnnotations$Name[which(annot$Type ==
        "I")]
    infiniumII <- probeAnnotations$Name[which(annot$Type ==
        "II")]
    cat("\t Quantile normalization of samples: separated and 'robust' quantile normalization for Infinium probes I and II through probe categories (reference quantiles computed from filtered Infinium I probes only and for different categories of probe annotations). \n")
    if (!is.null(annotCategory)) {
        if (annotCategory == "relationToCpG") {
            index <- which(is.element(colnames(annot),
                c("Name", "Relation_to_Island")))
            annot <- annot[, index]
        }
        if (annotCategory == "relationToSequence") {
            index <- which(is.element(colnames(annot),
                c("Name", "UCSC_RefGene_Group")))
            annot <- annot[, index]
        }
        if (!is.element(annotCategory, c("relationToCpG",
            "relationToSequence"))) {
            print("WARNINGS: probe annotation category must be one of 'relationToCpG' or 'relationToSequence'.")
            return("WARNINGS: probe annotation category must be on of 'relationToCpG' or 'relationToSequence'.")
        }
    }
    else {
        print("WARNING ! You have to specify an annotation type for probe categories based normalization ('relationToCpG' or 'relationToSequence').")
        return("WARNING ! You have to specify an annotation type for probe categories based normalization ('relationToCpG' or 'relationToSequence').")
    }
    print("ok")
    data.norm <- robustQuantileNorm_Illumina450K(data = beta,
        infiniumI = infiniumI, infiniumII = infiniumII, detect.pval = detect.pval,
        detect.pval.threshold = quantile.norm.pvalThreshold,
        annotations = annot)
    names(data.norm) <- c("beta", "detection.pvalue")
    rm(beta, detect.pval)
    rm(infiniumI, infiniumII, annot)
    cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1],
        "x", dim(data.norm$beta)[2], "\n")
    cat("Dimension of normalized detection p-values matrix: ",
        dim(data.norm$detection.pvalue)[1], "x", dim(data.norm$detection.pvalue)[2],
        "\n")
    return(data.norm)
}


robustQuantileNorm_Illumina450K.probeCategories
function (data.infiniumI, data.infiniumII, detect.pval.infiniumI,
    annotations, threshold, verbose = TRUE)
{
    res.norm <- list()
    data.detectPvalRemoved.infiniumI <- dataDetectPval2NA(data.infiniumI,
        detect.pval.infiniumI, threshold)
    category <- uniqueAnnotationCategory(annotations[, 2])
    if (verbose)
        cat("\tProbe categories values: ", paste(category, sep = "",
            collapse = ", "), sep = "")
    for (i in 1:length(category)) {
        if (verbose)
            cat("\n\tFor ", category[i], ":\n")
        probes_i <- findAnnotationProbes(annotationValue = category[i],
            annotations, uniqueAnnot = TRUE)
        robustProbes_i <- findAnnotationProbes(annotationValue = category[i],
            annotations, uniqueAnnot = TRUE)
        if (verbose)
            cat("\t\tNb probes for category", category[i], ":",
                length(probes_i), "\n")
        probeIndexI_i <- which(is.element(rownames(data.infiniumI),
            probes_i))
        if (verbose)
            cat("\t\tNb Inf_I probes for category ", category[i],
                ": ", length(probeIndexI_i), "\n")
        probeIndexII_i <- which(is.element(rownames(data.infiniumII),
            probes_i))
        if (verbose)
            cat("\t\tNb Inf_II probes for category ", category[i],
                ": ", length(probeIndexII_i), "\n")
        probe.detectPvalRemoved.IndexI_i <- which(is.element(rownames(data.detectPvalRemoved.infiniumI),
            robustProbes_i))
        if (verbose)
            cat("\t\tNb filtred Inf_I' probes for category ",
                category[i], ": ", length(probe.detectPvalRemoved.IndexI_i),
                "\n")
        data.infiniumI_i <- data.infiniumI[probeIndexI_i, ]
        data.infiniumII_i <- data.infiniumII[probeIndexII_i,
            ]
        data.detectPvalRemoved.infiniumI_i <- data.detectPvalRemoved.infiniumI[probe.detectPvalRemoved.IndexI_i,
            ]
        probeID.infiniumI_i <- rownames(data.infiniumI_i)
        probeID.infiniumII_i <- rownames(data.infiniumII_i)
        ref.quantiles.detectPvalRemoved.infI_i <- referenceQuantiles(data.detectPvalRemoved.infiniumI_i)
        indexRefQuantilesNA <- which(is.na(ref.quantiles.detectPvalRemoved.infI_i))
        if (length(indexRefQuantilesNA) > 0 && length(indexRefQuantilesNA) <
            length(ref.quantiles.detectPvalRemoved.infI_i)) {
            ref.quantiles.detectPvalRemoved.infI_i <- ref.quantiles.detectPvalRemoved.infI_i[-indexRefQuantilesNA]
        }
        if (length(indexRefQuantilesNA) == length(ref.quantiles.detectPvalRemoved.infI_i)) {
            print(paste("WARNING ! nb of 'NA' in ref.quantiles from dataI' for ",
                category[i], "equals the nb of quantiles (",
                length(indexRefQuantilesNA), ") !!!!", sep = ""))
            return(paste("WARNING ! nb of 'NA' in ref.quantiles from dataI' for ",
                category[i], "equals the nb of quantiles (",
                length(indexRefQuantilesNA), ") !!!!", sep = ""))
        }
        if ((length(indexRefQuantilesNA)/length(ref.quantiles.detectPvalRemoved.infI_i)) >
            0.8) {
            print(paste("WARNING !! nb of 'NA' in ref.quantiles from dataI' for ",
                category[i], "represents more than 80% of the ref.quantiles (",
                (length(indexRefQuantilesNA)/length(ref.quantiles.detectPvalRemoved.infI_i) *
                  100), ") !!!!", sep = ""))
        }
        ref.quantiles.infI_i <- adaptRefQuantiles(ref.quantiles.detectPvalRemoved.infI_i,
            dim(data.infiniumI_i)[1])
        ref.quantiles.infII_i <- adaptRefQuantiles(ref.quantiles.detectPvalRemoved.infI_i,
            dim(data.infiniumII_i)[1])
        data.infiniumI.norm_i <- normalize.quantiles2(data.infiniumI_i,
            ref.quantiles.infI_i)
        rownames(data.infiniumI.norm_i) <- probeID.infiniumI_i
        data.infiniumII.norm_i <- normalize.quantiles2(data.infiniumII_i,
            ref.quantiles.infII_i)
        rownames(data.infiniumII.norm_i) <- probeID.infiniumII_i
        res.norm$data.infiniumI.norm <- rbind(res.norm$data.infiniumI.norm,
            data.infiniumI.norm_i)
        res.norm$data.infiniumII.norm <- rbind(res.norm$data.infiniumII.norm,
            data.infiniumII.norm_i)
    }
    return(rbind(res.norm$data.infiniumI.norm, res.norm$data.infiniumII.norm))
}
