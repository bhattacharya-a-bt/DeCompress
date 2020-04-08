findInformSet <- function(yref,
                          method,
                          n_genes){

    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }


    if (method == 'linearity'){
        yref_need = linCor(yref)
    }

    if (method == 'variance'){
        yref_need = vardecomp(yref,
                              n_genes)
    }

    return(yref_need)


}
