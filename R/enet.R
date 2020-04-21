#' Elastic-net for compressed sensing
#'
#' The function runs elastic net to train compression for a single gene
#'
#' @param need vector, numeric vector of gene that is needed
#' @param train matrix, numeric expression matrix of training data
#' @param seed numeric, random seed
#' @param alpha numeric, between 0 and 1
#'
#' @return list with coefficients and predictive R2
#'
#' @importFrom glmnet cv.glmnet
#'
#' @export
enet <- function(need,
                 train,
                 seed = 1218,
                 alpha){

    if (alpha < 0 | alpha > 1){
        stop('alpha is not between 0 and 1')
    }


    if (class(need) != 'numeric'){
        stop('need is not a numeric vector')
    }


    if (class(train) != 'matrix'){
        stop('train is not a numeric matrix')
    }

    set.seed(seed)
    lar = glmnet::cv.glmnet(x = t(train),
                            y = need,
                            alpha = alpha,
                            k = 5)

    return(list(coef = as.numeric(coef(lar,s='lambda.min')[-1,1]),
                r2 = lar$glmnet.fit$dev.ratio[which(lar$glmnet.fit$lambda ==
                                                        lar$lambda.min)]))

}
