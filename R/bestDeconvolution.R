#' Use NMF repeatedly with TOAST framework to find cell-type gene signatures
#' and do deconvolution
#'
#' The function finds cell-type gene signatures and does deconvolution
#' using the TOAST framework with NMF
#'
#' @param yref matrix, numeric expression matrix
#' @param iters numeric, number of interations
#' @param pval numeric, p-value cutoff
#' @param n.types integer, number of cell-types
#' @param scree character, method to estimate n.types if n.types is NULL
#' @param log logical, T/F if yref is in log-scale
#'
#' @return list with cell-type proportions and expression profiles
#'
#' @export
bestDeconvolution <- function(yref,
                              n.types = NULL,
                              scree = c('drop','cumvar','residual'),
                              logTransform = F,
                              known.props = NULL){


    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    if (logTransform){
        yref = log2(yref+1)
    }

    row.means = rowSums(yref)
    yref = yref[row.means > 0,]

    if (is.null(n.types)){
        n.types = findNumberCells(yref,scree = scree)
    }

    geneList = rownames(yref)

    ### TOAST + NMF
    toast.nmf <- csDeCompress(Y_raw = yref,
                              K = n.types,
                              nMarker = min(1000,nrow(yref)),
                              FUN = nmfOut,
                              TotalIter = 30)

    require(NMF)
    fin.nmf = nmf(x = yref[toast.nmf$finalsigs,],
                  rank = n.types)
    ppp = t(coef(fin.nmf))
    ppp = t(apply(ppp,1,function(c) c/sum(c)))
    nmf.res = list(prop = ppp,
                   sig = basis(fin.nmf))


    ### LINSEED
    linseed.rs = linCor(yref = yref,
                        iters = 100,
                        pval = .01,
                        n.types = n.types,
                        scree = 'drop',
                        logTransform = F)
    names(linseed.rs) = names(nmf.res)
    linseed.rs$prop = t(linseed.rs$prop)

    ### CellDistinguisher
    cd <- CellDistinguisher::gecd_CellDistinguisher(yref,
                                 genesymb=geneList,
                                 numCellClasses=n.types,
                                 minDistinguisherAlternatives=1,
                                 maxDistinguisherAlternatives=100,
                                 minAlternativesLengthsNormalized=0.5,
                                 expressionQuantileForFilter=0.999,
                                 expressionConcentrationRatio=0.333,
                                 verbose=0)
    CellDist.deconv <-
        tryCatch(CellDistinguisher::gecd_DeconvolutionByDistinguishers(
        as.matrix(yref),
        cd$bestDistinguishers,
        nonNegativeOnly = FALSE,
        convexSolution = FALSE,
        verbose = 0),
        error = function(e) return(list(sampleCompositions =
                                            matrix(rep(1/n.types,
                                                       n.types*ncol(yref)),
                                                   ncol=n.types))))
    if (length(CellDist.deconv) > 1){
        CellDist.rs = list(prop = t(CellDist.deconv$sampleCompositions),
                           sig = CellDist.deconv$cellSubclassSignatures)
    } else {CellDist.rs = nmf.res}

    ### DeconICA
    deconica.deconv = deconica::run_fastica(yref,
                                            overdecompose = F,
                                            with.names = F,
                                            gene.names = geneList,
                                            samples = colnames(yref),
                                            n.comp = n.types,
                                            R = T)
    deconica.rs = list(prop = t(deconica.deconv$A),
                       sig = deconica.deconv$S)

    all.res = list(nmf.res,
                   linseed.rs,
                   CellDist.rs,
                   deconica.rs)
    errors = vector(mode = 'numeric',
                    length = length(all.res))
    for (i in 1:length(errors)){
        errors[i] = calculateError(all.res[[i]]$prop,
                                   all.res[[i]]$sig,
                                   trueProp = known.props,
                                   yref = yref)
    }




    return(all.res[[which.min(errors)]])
}
