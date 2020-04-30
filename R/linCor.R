#' Use mutual linearity to find cell-type gene signatures and do deconvolution
#'
#' The function finds cell-type gene signatures and does deconvolution
#' using the linseed pacakge
#'
#' @param yref matrix, numeric expression matrix
#' @param iters numeric, number of interations
#' @param pval numeric, p-value cutoff
#' @param n.types integer, number of cell-types
#' @param scree character, method to estimate n.types if n.types is NULL
#'
#' @return list with coefficients and predictive R2
#'
#' @export
linCor <- function(yref,
                   iters = 1000,
                   pval = .01,
                   n.types = NULL,
                   scree = c('drop','cumvar','residual')){

    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    lo <- linseed::LinseedObject$new(yref)
    lo$calculatePairwiseLinearity()
    lo$calculateSpearmanCorrelation()
    lo$calculateSignificanceLevel(iters)
    lo$filterDatasetByPval(pval)

    if (is.null(n.types)){
        s = svd(lo$exp$filtered$norm)
        if (scree == 'residual'){
        p = 10
        y <- lapply(1:(p-1),function(i)
            try(structureMeasure(s,lo$exp$filtered$norm,i)))
        y2 <- unlist(y[sapply(y,
                              function(x) !inherits(x, "try-error"))])
        y2.ind <- which(sapply(y,
                               function(x) !inherits(x, "try-error")))
        n.types = y2.ind[which.min(y2)]
        }
        if (scree == 'cumvar'){
            dii = s$d^2
            cumvar = cumsum(dii/sum(dii))
            n.types = min(which(cumvar >= .8))
        }
        if (scree == 'drop'){
            n.types = kaiser(s$d^2)
        }
    }

    lo$setCellTypeNumber(n.types)
    lo$project("full")
    lo$project("filtered")
    lo$smartSearchCorners(dataset="filtered", error="norm")
    lo$deconvolveByEndpoints()


    return(list(prop = lo$proportions,
                sig = lo$signatures))

}
