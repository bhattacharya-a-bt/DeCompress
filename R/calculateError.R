calculateError <- function(estProp,
                           estSig,
                           trueProp = NULL,
                           yref){

    if (!is.null(trueProp)){
        return(mseError(estProp,trueProp,rand=F))
    } else {
            yref = yref[rownames(estSig),]
            return(mseError(yref,estSig %*% t(estProp),rand=F))
        }

}
