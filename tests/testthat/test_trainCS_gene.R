test_that("compression training works", {

    dir = system.file('extdata',package = 'DeCompress')
    gtex = readRDS(paste0(dir,'/example_data.RDS'))
    ref = gtex$true.mixed.expression
    target = gtex$observed.mixed

    compress = trainCS_gene(need = ref[1,],
                            train = as.matrix(ref[rownames(target)[1:200],]),
                            seed = 1218,
                            method = c('lar',
                                       'lasso',
                                       'l1'),
                            n.cores = 1,
                            lambda = .1)

    expect_equal(length(compress),2)
    expect_equal(length(compress$coef),200)


})
