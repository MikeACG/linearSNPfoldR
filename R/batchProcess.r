#' @export
SNPfoldR <- function(trId, trIdx, rna, trsLength, bppsDir, bppsSuffix) {

    # infer the path to this RNA results based on the length of all rnas, and the index, name and seq of query
    pRNAchanges <- getRNAchanges(rna)
    changesBatchDirs <- trIndex2batchDirs(trIdx, trsLength)
    bppsFileWt <- paste0(bppsDir, changesBatchDirs[1], trId, "_WT", bppsSuffix)
    bppsFilesMuts <- paste0(bppsDir, changesBatchDirs[-1], trId, "_", pRNAchanges, bppsSuffix)

    # load wild type bpps
    n <- length(pRNAchanges) / 3
    bppsWt <- file2bppMat(bppsFileWt, n)

    # compute PCC for each possible RNA change against wild type
    pccs <- sapply(bppsFilesMuts, function(f) halvorsenPCC(bppsWt, file2bppMat(f, n)))

    # organize results, they are already sorted by position of change
    rnaChangePos <- as.integer(gsub("^[A-Z]|[A-Z]$", "", pRNAchanges))
    changeTo <- gsub(".*[0-9]", "", pRNAchanges)
    SNPfoldDt <- data.table::data.table(RNA_Change = pRNAchanges, PCC = pccs, Change_Index = rnaChangePos, Change_To = changeTo)
    #SNPfoldDt <- SNPfoldDt[order(SNPfoldDt$Change_Index, SNPfoldDt$Change_To), ]

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
    t2 <- t1 * ((3 * rnaLens) + 1)

    # output total time in hours
    return(sum(t2) / 60)

}

trIndex2batchDirs <- function(trIndex, trsLength) {

    maxFilesPerDir <- 25000

    # get the amount of results that have to be saved for this transcript
    nresults <- (trsLength[trIndex] * 3) + 1

    # compute the amount of results that were stored before results for this transcript
    accResults <- 0
    if (trIndex > 1) {

        beforeLengths <- trsLength[1:(trIndex - 1)]
        accResults <- sum((beforeLengths * 3) + 1)

    }

    # get the corresponding batch folder for each of the results of this transcript
    resultsIdxs <- (accResults + 1):(accResults + nresults)
    batchDirs <- paste0("batch", ceiling(resultsIdxs / maxFilesPerDir), "/")

    return(batchDirs)

}

#' @export
makeAll2foldFile <- function(trsRNA, trsId, all2foldFilePath) {

    if (file.exists(all2foldFilePath)) unlink(all2foldFilePath)
    file.create(all2foldFilePath)

    # drop missing values
    isValid <- !is.na(trsRNA)
    trsRNA <- trsRNA[isValid]
    trsId <- trsId[isValid]

    trsLength <- sapply(trsRNA, nchar)
    ntrs <- length(trsRNA)
    for (i in 1:ntrs) {

        if (i %% 100 == 0) cat(paste0(i, "/", ntrs, "...\n"))

        pRNAchanges <- c("WT", getRNAchanges(trsRNA[i]))
        pMutRNAs <- c(trsRNA[i], getMutRNAs(trsRNA[i], pRNAchanges[-1]))
        foldDt <- data.table::data.table(mut = pRNAchanges, rna = pMutRNAs)

        foldDt$id <- trsId[i]
        foldDt$batchdir <- trIndex2batchDirs(i, trsLength)
        foldDt$filename <- paste0(foldDt$id, "_", pRNAchanges)
        foldDt$path <- paste0(foldDt$batchdir, foldDt$filename)

        data.table::fwrite(foldDt, all2foldFilePath, sep = "\t", col.names = FALSE, 
            quote = FALSE, append = TRUE, na = "NA"
        )

    }

}

