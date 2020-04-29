nmfOut <- function(Y,
                   K){

    outY = NMF::nmf(Y,rank = K)
    return(t(coef(outY)))

}
