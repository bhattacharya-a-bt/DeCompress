expandTarget <- function(target,
                         compression_mat){

    if (class(target) != 'matrix'){
        stop("matrix not supplied in target")
    }

    if (class(compression_mat) != 'matrix'){
        stop("matrix not supplied in compression_mat")
    }

    if (nrow(compression) == nrow(target)){
        return(t(compression) %*% target)
    } else {return(compression %*% target)}

}
