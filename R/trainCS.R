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
    compression = pbapply::pbapply(yref_need,
                                   MARGIN = 1,
                                   trainCS_gene,
                                   train = yref,
                                   seed=seed,
                                   method = method)
    }

    if (par){
        yref_need_list = split(yref_need,
                               rep(1:nrow(yref_need),
                                   each = ncol(yref_need)))
        compression = parallel::mclapply(yref_need_list,
                                         trainCS_gene,
                                         train = yref,
                                         seed=seed,
                                         method = method,
                                         mc.cores = n.cores)
    }

    compression_mat = sapply(compression,function(x) x$coef)
    r2Vec = sapply(compression,function(x) x$r2)

    return(list(compression.matrix = compression_mat,
                r2 = r2Vec))

}
