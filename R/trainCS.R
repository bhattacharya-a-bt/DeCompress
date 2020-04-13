trainCS <- function(yref,
                    yref_need,
                    target,
                    seed = 1218,
                    method = c('lar',
                               'lasso',
                               'enet',
                               'ridge',
                               'l1',
                               'TV',
                               'l2')){

    if (class(yref) != 'matrix'){
        stop("matrix not supplied in yref")
    }

    if (class(yref_need) != 'matrix'){
        stop("matrix not supplied in yref_need")
    }

    if (class(target) != 'matrix'){
        stop("matrix not supplied in target")
    }

    compression = pbapply::pbapply(yref_need,
                                   MARGIN = 1,
                                   trainCS_gene,
                                   train = yref,
                                   seed=seed,
                                   method = method)

    return(compression)

}
