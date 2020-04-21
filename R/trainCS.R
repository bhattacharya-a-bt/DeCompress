trainCS <- function(yref,
                    yref_need,
                    seed = 1218,
                    method = c('lar',
                               'lasso',
                               'enet',
                               'ridge',
                               'l1',
                               'TV',
                               'l2'),
                    par = T,
                    n.cores){

    if (class(yref) != 'matrix'){
        stop("matrix not supplied in yref")
    }

    if (class(yref_need) != 'matrix'){
        stop("matrix not supplied in yref_need")
    }

    if (!par){
    compression = apply(yref_need,
                        MARGIN = 1,
                        trainCS_gene,
                        train = yref,
                        seed=seed,
                        method = method)
    }

    if (par){
        future::plan(multiprocess,workers = n.cores)
        compression = future.apply::future_apply(yref_need,
                                                 MARGIN = 1,
                                                 trainCS_gene,
                                                 train = yref,
                                                 method = method,
                                                 seed = seed)
    }

    compression_mat = sapply(compression,function(x) x$coef)
    r2Vec = sapply(compression,function(x) x$r2)

    return(list(compression.matrix = compression_mat,
                r2 = r2Vec))

}
