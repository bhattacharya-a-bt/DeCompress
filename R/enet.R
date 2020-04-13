enet <- function(need,
                 train,
                 seed = 1218){

    lar = tryCatch(glmnet::cv.glmnet(x = t(train),
                                     y = need,
                                     alpha = .5,
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
