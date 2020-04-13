enet <- function(need,
                 train,
                 seed = 1218,
                 alpha){

    set.seed(seed)
    lar = glmnet::cv.glmnet(x = t(train),
                            y = need,
                            alpha = alpha,
                            k = 5)

    return(list(coef = as.numeric(coef(lar,s='lambda.min')[-1,1]),
                r2 = lar$glmnet.fit$dev.ratio[which(lar$glmnet.fit$lambda ==
                                                        lar$lambda.min)]))

}
