#' Expand targetted panel to larger feature space
#'
#' The function multiplies the targetted panel with the compression
#' matrix for the larger feature space
#'
#' @param target matrix, numeric expression matrix from targetted panel
#' @param compression_mat matrix, compression matrix
#'
#' @return expanded expression matrix
#'
#' @export
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
