csDeCompress <- function(Y_raw,
                         K,
                         FUN = nmfOut,
                         nMarker = 1000,
                         InitMarker = NULL,
                         TotalIter = 30,
                         bound_negative = FALSE) {

    if (is(Y_raw, "SummarizedExperiment")) {
        se <- Y_raw
        Y_raw <- assays(se)$counts
    } else if (!is(Y_raw, "matrix")) {
        stop("Y_raw should be a matrix or a SummarizedExperiment object!")
    }

    if (is.null(rownames(Y_raw))) {
        row.names(Y_raw) <- seq(nrow(Y_raw))
    }
    if (is.null(InitMarker)) {
        if (nrow(Y_raw) < 2000) {
            InitMarker <- TOAST::findRefinx(Y_raw, nmarker = nMarker)
        } else {
            tmp <- TOAST::findRefinx(Y_raw, nmarker = nMarker*2)
            InitMarker <- tmp[nMarker+1:nMarker]
        }
    } else {
        if (sum(!(InitMarker %in% rownames(Y_raw))) > 0) {
            stop("Discrepancy between
                    InitMarker and the row names of Y_raw!")
        }
    }

    allProp <- list()
    allRMSE <- rep(0, TotalIter + 1)

    Y <- Y_raw[InitMarker, ]
    Prop0 <- FUN(Y, K)
    allProp[[1]] <- Prop0

    out_all <- csSAM::csfit(Prop0, t(Y_raw))
    prof <- t(out_all$ghat)
    tmpmat <- prof %*% t(Prop0)
    allRMSE[1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))

    message("+========================================+")
    message("+======= Total iterations = ",
            TotalIter, " ==========+")

    for (i in seq_len(TotalIter)) {
        message("Current iter = ", i)

        updatedInx <- TOAST::DEVarSelect(Y_raw, Prop0,
                                         nMarker, bound_negative)
        Y <- Y_raw[updatedInx, ]
        Prop0 <- FUN(Y, K)
        allProp[[i + 1]] <- Prop0

        out_all <- csSAM::csfit(Prop0, t(Y_raw))
        prof <- t(out_all$ghat)
        tmpmat <- prof %*% t(Prop0)
        allRMSE[i + 1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))
    }

    min_idx <- which.min(allRMSE)
    Prop0 <- allProp[[min_idx]]

    return(list(allRMSE = allRMSE,
                allProp = allProp,
                estProp = Prop0,
                finalsigs = updatedInx
                )
           )
}