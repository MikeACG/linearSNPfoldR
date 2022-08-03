#' @export
SNPfoldR <- function(rnaDir, bppFilesName) {

    # get paths to directories with mutated foldings
    rnaDirFiles <- list.files(rnaDir)
    rnaMutBatchDirs <- paste0(rnaDir, rnaDirFiles[grep("batch", rnaDirFiles)], "/")
    rnaMutDirs <- lapply(rnaMutBatchDirs, list.files)
    rnaMutPaths <- unlist(mapply(function(b, m) paste0(b, m, "/", bppFilesName), 
        rnaMutBatchDirs, rnaMutDirs, SIMPLIFY = FALSE
    ))

    # get sequence length from the mutations information
    rnaMuts <- unlist(rnaMutDirs)
    rnaChangePos <- as.integer(gsub("^[A-Z]|[A-Z]$", "", rnaMuts))
    n <- max(rnaChangePos)

    # load wild type bpps
    bppFileWt <- paste0(rnaDir, "WT/", bppFilesName)
    bppWt <- file2bppMat(bppFileWt, n)

    # compute PCC for each possible RNA change against wild type
    pccs <- sapply(rnaMutPaths, function(f) halvorsenPCC(bppWt, file2bppMat(f, n)))

    # organize results sorted by position of change and A, C, G, U order
    changeTo <- gsub(".*[0-9]", "", rnaMuts)
    SNPfoldDt <- data.table::data.table(RNA_Change = rnaMuts, PCC = pccs, Change_Index = rnaChangePos, Change_To = changeTo)
    SNPfoldDt <- SNPfoldDt[order(SNPfoldDt$Change_Index, SNPfoldDt$Change_To), ]

    return(SNPfoldDt)

}

makeDecompressCmd <- function(fileName) {

    if (grepl(".xz$", fileName)) return(paste("xz -dc", fileName))

    return("notCompressed")

}

#' @export
computeTime <- function(rnas) {

    # parameters according to Huang et al on a Linux machine with 2.90 GHz Intel i9-7920X CPU
    m <- 1.3
    l <- 32753

    rnaLens <- sapply(rnas, nchar)

    # time taken to fold each sequence once
    t1 <- (rnaLens * m) / l

    # each sequence must be folded 3n + 1 times
    t2 <- t1 * ((3 * rnaLens) + rnaLens)

    # output total time in hours
    return(sum(t2) / 60)

}

ensembl2folder <- function(id, granularity) {

    ensPrefix <- gsub("[0-9]*", "", id)
    nid <- as.numeric(gsub("ENS.", "", id)) # trim letters and leading zeros
    rnid <- format(nid - (nid %% granularity), scientific = FALSE) # floor to nearest granularity multiple
    n <- nchar(rnid) # n of numbers in the id (without leading zeros)
    pad <- paste(rep("0", pmax(11, n) - n), collapse = '') # number of leading zeros, have to complete eleven numbers at least
    ll <- paste0(ensPrefix, pad, rnid) # lower limit (inclusive) of folder
    ul <- paste0(ensPrefix, pad, as.numeric(rnid) + (granularity - 1)) # upper limit of folder (non-inclusive)
    
    return(paste0(ll, "_", ul, "/"))

}

mut2folder <- function(muts, granularity) {

    mutsSplit <- split(muts, ceiling(seq_along(muts) / granularity))
    dirs <- rep(paste0("batch", 1:length(mutsSplit), "/"), sapply(mutsSplit, length))

    return(dirs)

}

#' @export
makeAll2foldFile <- function(trsRNA, trsEnsembl, all2foldFilePath) {

    maxTrsPerFolder <- 1000
    maxMutsPerFolder <- 500
    if (file.exists(all2foldFilePath)) unlink(all2foldFilePath)
    file.create(all2foldFilePath)

    # drop missing values
    isValid <- !is.na(trsRNA)
    trsRNA <- trsRNA[isValid]
    trsEnsembl <- trsEnsembl[isValid]

    ntrs <- length(trsRNA)
    for (i in 1:ntrs) {

        if (i %% 100 == 0) cat(paste0(i, "/", ntrs, "...\n"))

        pRNAchanges <- getRNAchanges(trsRNA[i])
        pMutRNAs <- getMutRNAs(trsRNA[i], pRNAchanges)
        foldDt <- data.table::data.table(mut = pRNAchanges, rna = pMutRNAs)

        foldDt$ensembl <- trsEnsembl[i]
        foldDt$trdir <- ensembl2folder(trsEnsembl[i], maxTrsPerFolder)
        foldDt$mutdir <- mut2folder(foldDt$mut, maxMutsPerFolder)
        foldDt$path <- paste0(foldDt$trdir, foldDt$ensembl, "/", foldDt$mutdir, foldDt$mut, "/")

        # add wild type in its own folder
        foldDt <- rbind(data.table::data.table(mut = "WT", rna = trsRNA[i], 
            ensembl = trsEnsembl[i], trdir = foldDt$trdir[1], mutdir = "WT/",
            path = paste0(foldDt$trdir[1], trsEnsembl[i], "/", "WT/")
        ), foldDt)

        data.table::fwrite(foldDt, all2foldFilePath, sep = "\t", col.names = FALSE, 
            quote = FALSE, append = TRUE
        )

    }

}
