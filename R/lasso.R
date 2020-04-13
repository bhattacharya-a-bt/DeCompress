lasso <- function(need,
                  train,
                  seed = 1218){

    lar = tryCatch(lars::cv.lars(x = t(train),
                                 y = need,
                                 type = 'lasso',
                                 plot.it = F,
                                 K = 5),
                   error = function(e) {return(NULL)})
    if (!is.null(lar)){
        lar.reg = lars(x = t(train),
                       y = need,
                       type = 'lasso')
        lar.coef = as.numeric(coef(lar.reg,
                                   s = lar$index[which.min(lar$cv)],
                                   mode='fraction'))
    } else {
        lar.coef = NULL
        lar = list(cv = var(need))
    }

    return(list(coef = lar.coef,
                r2 = 1 - min(lar$cv)/var(need)))

}
