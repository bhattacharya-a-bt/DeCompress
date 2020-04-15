parL1 <- function(phi,
                  y,
                  T,
                  x0,
                  lambda = 0.1,
                  n.cores){

    cl <- parallel::makeCluster(n.cores)
    setDefaultCluster(cl)
    qqq = optimParallel::optimParallel(par = x0,
                                       fn = R1magic::objectiveL1,
                                       y = y,
                                       T = T,
                                       lambda = lambda)
    setDefaultCluster(cl=NULL)
    stopCluster(cl)

    return(qqq)

}
