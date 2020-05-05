nmfOut <- function(Y,
                   K){

    require(NMF)
    outY = nmf(Y,rank = K)
    return(t(coef(outY)))

}
