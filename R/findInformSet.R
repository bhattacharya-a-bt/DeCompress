#' Select the number of cell types using SVD methods
#'
#' The function estimates the number of cell-types using one
#' of three SVD methods on the input expression matrix
#'
#' @param yref matrix, numeric expression matrix
#' @param method character, variance/linearity for TOAST or LINSEED
#' @param n_genes integer, number of CTS genes needed
#' @param n.types integer, number of compartments
#' @param scree character, method to estimate n.types
#'
#' @return numeric matrix of n_genes CTS genes
#'
#' @export
findInformSet <- function(yref,
                          method = c('variance','linearity'),
                          n_genes = 1000,
                          n.types = NULL,
                          scree = 'cumvar'){


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
