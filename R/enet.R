enet <- function(need,
                 train,
                 seed = 1218,
                 alpha){

    set.seed(seed)
    lar = tryCatch(glmnet::cv.glmnet(x = t(train),
                                     y = need,
                                     alpha = alpha,
                                     k = 5),
                   error = function(e) {return(NULL)})
    if (is.null(lar)){
        lar.coef = NULL
        lar = list(cv = var(need))
        }

    return(list(coef = lar.coef,
                r2 = lar$glmnet.fit$dev.ratio[which(lar$glmnet.fit$lambda ==
                                                        lar$lambda.min)]))

}
