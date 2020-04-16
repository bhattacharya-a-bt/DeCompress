parMagic <- function(train,
                     need,
                     lambda = 0.1,
                     seed = 1218,
                     obj,
                     n.cores){

    require(R1magic)
    set.seed(seed)
    Tt = diag(nrow(train))
    p = matrix(0,nrow(train),1)

    if (obj == 'l1'){
        func = objectiveL1
    }
    if (obj == 'l2'){
        func = objectiveL2
    }
    if (obj == 'TV'){
        func = objective1TV
    }


    cl <- parallel::makeCluster(n.cores)
    parallel::setDefaultCluster(cl)
    qqq = optimParallel::optimParallel(par = as.vector(p),
                                       fn = func,
                                       phi = t(train),
                                       y = need,
                                       T = Tt,
                                       lambda = lambda,
                                       parallel=list(forward=F,
                                                     loginfo=F),
                                       verbose = T,
                                       method = 'L-BFGS-B')

    ### DO CV
    pred = vector('numeric',
                  length(need))
    train.folds = caret::createFolds(y = need,
                                     k=5,
                                     returnTrain = T)

    for (i in 1:length(train.folds)){

        new.train = train[,train.folds[[i]]]
        Tt = diag(nrow(new.train))
        p = matrix(0,nrow(new.train),1)
        new.mod = optimParallel::optimParallel(par = as.vector(p),
                                               fn = objectiveL1,
                                               phi = t(new.train),
                                               y = need[train.folds[[i]]],
                                               T = Tt,
                                               lambda = lambda,
                                               parallel=list(forward=F,
                                                             loginfo=F),
                                               verbose = T,
                                               method = 'L-BFGS-B')
        pred[-train.folds[[i]]] = new.mod$par %*% train[,-train.folds[[i]]]


    }

    parallel::setDefaultCluster(cl=NULL)
    parallel::stopCluster(cl)


    return(list(coef = qqq$par,
                r2 = cor(need,pred)^2))

}
