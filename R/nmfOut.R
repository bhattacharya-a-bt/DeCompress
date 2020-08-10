#' NMF function for TOAST
#'
#' The function runs NMF on an expression matrix with rank K
#'
#' @param Y matrix, numeric expression matrix
#' @param K integer, number of cell-types or rank of factorization
#'
#' @return coefficient matrix from NMF
#'
#'
#' @importFrom NMF nmf
#'
#' @export
nmfOut <- function(Y,
                   K){

    outY = NMF::nmf(Y,rank = K)
    return(t(coef(outY)))

}
