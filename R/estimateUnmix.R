estimateUnmix <- function(yref,
                          sigs,
                          shift.range = 10^(1:10),
                          logTransform = F){


    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    if (logTransform){
        yref = log2(yref+1)
    }

    row.means = rowSums(yref)
    yref = yref[row.means > 0,]

    geneList = rownames(yref)

    regBetweenRows <- function(shift,mat){
        expressions <- log(mat+shift+1)
        means <- rowMeans(expressions)
        sds <- apply(expressions,1,sd)
        return(as.numeric(abs(coef(lm(sds~means))[2])))
    }

    slope = sapply(shift.range,
                   regBetweenRows,
                   mat = yref)
    shift = shift.range[which.min(slope)]
    props = DESeq2::unmix(as.matrix(yref),
                          pure=as.matrix(sigs),
                          shift = shift)
    props = t(apply(props,1,function(v) v/sum(v)))

    return(list(prop = props,
                sig = sigs))
}
