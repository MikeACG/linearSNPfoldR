
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
    bpps <- setNames(data.table::fread(bppFile, colClasses = c("integer", "integer", "numeric")), c("nuc1", "nuc2", "bpp"))

    for (j in 1:nrow(bpps)) {

        nuc1 <- bpps$nuc1[j]
        nuc2 <- bpps$nuc2[j]
        bppMat[nuc1, nuc2] <- bpps$bpp[j]

    }

    return(bppMat)

}

#' @export
halvorsenPCC <- function(X, Y) {

    x <- apply(X, 2, sum)
    y <- apply(Y, 2, sum)

    bothNonzero <- x > 0 & y > 0
    pcc <- cor(x[bothNonzero], y[bothNonzero], method = "pearson")

    return(pcc) 

}
