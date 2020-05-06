mseError <- function(mat1,mat2,rand=F){

    if (ncol(mat1) > nrow(mat1)){
        mat1 = t(mat1)
    }

    if (ncol(mat2) > nrow(mat2)){
        mat2 = t(mat2)
    }

    if (nrow(mat1) != nrow(mat2)){
        stop('mat1 and mat2 have different numbers of rows')
    }

    mse.row = function(i,mat1,mat2){
        return(hydroGOF::mse(mat1[i,],mat2[i,]))
    }

    mse.all = function(mat1,mat2){
        return(mean(sapply(1:nrow(mat1),
                           mse.row,
                           mat1=mat1,
                           mat2=mat2)))
    }

    if (rand == F){
        return(mse.all(mat1,mat2))
    }

    if (rand == T){

        p = gtools::permutations(n = ncol(mat1),
                         r = ncol(mat1),
                         v = 1:ncol(mat1),
                         repeats.allowed=F)
        error = 100
        for (j in 1:nrow(p)){
            error = min(error,mse.all(mat1,mat2[,p[j,]]))
        }
        return(error)

    }

}
