linCor <- function(yref,
                   iters = 1000,
                   pval = .01){

    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    lo <- linseed::LinseedObject$new(yref)
    lo$calculatePairwiseLinearity()
    lo$calculateSpearmanCorrelation()
    lo$calculateSignificanceLevel(iters)
    lo$filterDatasetByPval(pval)

    return(yref[rownames(yref) %in%
                    rownames(lo$exp$filtered$norm),])

}
