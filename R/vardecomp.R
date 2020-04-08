vardecomp <- function(yref,
                      n_genes = 1000){

    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    return(yref[TOAST::findRefinx(yref,
                                  nmarker = n_genes),])

}
