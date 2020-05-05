iterateSigs <- function(yref,
                        n.types = NULL,
                        n.markers = 1000,
                        n.iter = 30,
                        scree = c('drop','cumvar','residual'),
                        log = T){


    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    if (!log){
        yref = log2(yref+1)
    }

    row.means = rowSums(yref)
    yref = yref[row.means > 0,]

    if (is.null(n.types)){
        n.types = findNumberCells(yref,scree = scree)
        }

    geneList = rownames(yref)
    rownames(yref) = NULL
    outRF1 <- csDeCompress(Y_raw = yref,
                           K = n.types,
                           nMarker = n.markers,
                           FUN = nmfOut,
                           TotalIter = n.iter)
    return(list(estProp = outRF1$estProp,
                geneSig = geneList[outRF1$finalsigs]))
}
