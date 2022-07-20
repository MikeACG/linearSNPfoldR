
#' @export
linearSNPfoldR <- function(rna, seqid) {

    program <- "/usr/local/bin/LinearPartition-1.0/linearpartition -r"
    if (!file.exists(seqid)) dir.create(seqid)

    cat("Folding wild type sequence...\n")
    wtFile <- paste(seqid, "wt.txt", sep = "/")
    wtCmd <- paste("echo", rna, "|", program, wtFile)
    system(wtCmd)

    cat("Folding possible mutated sequences")
    pRNAchanges <- getRNAchanges(rna)
    pMutRNAs <- getMutRNAs(rna, pRNAchanges)
    npChanges <- length(pRNAchanges)
    mutFiles <- paste0(seqid, "/", 1:npChanges, ".txt")
    for (i in 1:npChanges) {

        cat(paste0(i, "/", npChanges, "...\n"))

        mutCmd <- paste("echo", pMutRNAs[i], "|", program, mutFiles[i])
        system(mutCmd)

    }

    cat("Computing correlation coefficients...\n")
    n <- nchar(rna)
    bppWt <- file2bppMat(wtFile, n)
    pccs <- sapply(mutFiles, function(f) halvorsenPCC(bppWt, file2bppMat(f, n)))

    outDt <- data.table(RNA_Change = pRNAchanges, PCC = pccs) 
    unlink(seqid, recursive = TRUE)
    return(outDt)

}

#' @export
getMutRNAs <- function(rna, pRNAchanges) {

    changeIdxs <- as.integer(gsub("^[A-Z]|[A-Z]$", "", pRNAchanges))
    changeNucs <- gsub("[0-9]*|^.", "", pRNAchanges)

    pMutRNAs <- rep(rna, length(pRNAchanges))
    substr(pMutRNAs, changeIdxs, changeIdxs) <- changeNucs

    return(pMutRNAs)

}

#' @export
getRNAchanges <- function(rna) {

    nucs <- c("A", "C", "G", "U")

    rnaSplit <- unlist(strsplit(rna, ""))
    possibleMuts <- lapply(rnaSplit, function(nuc) nucs[!(nucs %in% nuc)])
    changes <- mapply(paste0, rnaSplit, 1:length(rnaSplit), possibleMuts, SIMPLIFY = FALSE)

    return(unlist(changes))

}

#' @export
file2bppMat <- function(bppFile, n) {

    bppMat <- matrix(0, nrow = n, ncol = n)

    # read the base pairing probabilities file according to format
    cmd <- makeDecompressCmd(bppFile)
    if (cmd == "notCompressed") {

        bpps <- try(setNames(data.table::fread(bppFile, colClasses = c("integer", "integer", "numeric")), c("nuc1", "nuc2", "bpp")), silent = TRUE)

    } else {

        bpps <- try(setNames(data.table::fread(cmd = cmd, colClasses = c("integer", "integer", "numeric")), c("nuc1", "nuc2", "bpp")), silent = TRUE)

    }

    # if file reading fails (e.g. Folding has no non-zero probabilities, this happens for sure when seq length is 1) return zeros matrix
    if (class(bpps)[1] == "try-error") return(bppMat)

    # fill in the non-zero pairing probabilities
    for (j in 1:nrow(bpps)) {

        nuc1 <- bpps$nuc1[j]
        nuc2 <- bpps$nuc2[j]
        bppMat[nuc1, nuc2] <- bpps$bpp[j]

    }

    return(bppMat)

}

makeDecompressCmd <- function(fileName) {

    if (grepl(".xz$", fileName)) return(paste("xz -dc", fileName))

    return("notCompressed")

}

#' @export
halvorsenPCC <- function(X, Y) {

    x <- apply(X, 2, sum)
    y <- apply(Y, 2, sum)

    bothNonzero <- x > 0 & y > 0
    if (sum(bothNonzero) == 0) return(1)
    pcc <- cor(x[bothNonzero], y[bothNonzero], method = "pearson")

    return(pcc) 

}

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
