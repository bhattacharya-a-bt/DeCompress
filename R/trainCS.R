trainCS <- function(yref,
                    yref_need,
                    target,
                    method){

    if (class(yref) != 'matrix'){
        stop("matrix not supplied in yref")
    }

    if (class(yref_need) != 'matrix'){
        stop("matrix not supplied in yref_need")
    }

    if (class(target) != 'matrix'){
        stop("matrix not supplied in target")
    }



}
