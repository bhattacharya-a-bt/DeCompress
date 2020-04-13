lar <- function(need,
                train,
                seed = 1218){
    set.seed(seed)
    lar = tryCatch(lars::cv.lars(x = t(train),
                                 y = need,
                                 type = 'lar',
                                 plot.it = F,
                                 K = 5,
                                 use.Gram = T),
                   error = function(e) {return(NULL)})
    if (!is.null(lar)){
        lar.reg = lars(x = t(train),
                       y = need,
                       type = 'lar')
        lar.coef = as.numeric(coef(lar.reg,
                                   s = lar$index[which.min(lar$cv)],
                                   mode='step'))
        } else {
            lar.coef = NULL
            lar = list(cv = var(need))
            }

    return(list(coef = lar.coef,
                r2 = 1 - min(lar$cv)/var(need)))

}
