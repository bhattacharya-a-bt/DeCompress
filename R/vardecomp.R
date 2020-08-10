#' Select CTS genes with variance proporties
#'
#' The function runs TOAST's findRefinx function to find
#' a set of top cell-type specific genes
#'
#' @param yref matrix, numeric matrix of reference expression
#' @param n_genes integer, number of genes requested
#'
#' @return numeric matrix of n_genes CTS genes
#'
#' @importFrom TOAST findRefinx
#'
#' @export
vardecomp <- function(yref,
                      n_genes = 1000){

    if (class(yref) != c('matrix')){
        stop("matrix not supplied in yref")
    }

    return(yref[TOAST::findRefinx(yref,
                                  nmarker = n_genes),])

}
