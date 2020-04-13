trainCS_gene <- function(need,
                         train,
                         method){

    if (!class(need) %in% c('vector','numeric')){
        stop("provide a numeric vector for need")
    }

    if (class(train) != 'matrix'){
        stop("matrix not supplied for train")
    }



}
