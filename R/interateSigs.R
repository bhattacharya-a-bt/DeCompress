iterateSigs <- function(yref,
                        n.types = NULL,
                        n.markers = 1000,
                        n.iter = 30){


    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    if (is.null(n.types)){
        s = svd(yref)
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
