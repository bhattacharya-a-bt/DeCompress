test_that("deconvolution works", {

    dir = system.file('extdata',package = 'DeCompress')
    gtex = readRDS(paste0(dir,'/example_data.RDS'))
    ref = gtex$true.mixed.expression
    target = gtex$observed.mixed


    decompress.res = bestDeconvolution(target,
                                       n.types = 4,
                                       scree = 'cumvar',
                                       logTransform = F,
                                       known.props = NULL,
                                       methods = c('TOAST',
                                                   'linseed',
                                                   'celldistinguisher'))

    expect_equal(class(decompress.res),'list')
    expect_equal(length(decompress.res), 2)
    expect_equal(class(decompress.res[[1]]), c('matrix','array'))
    expect_equal(class(decompress.res[[2]]), c('matrix','array'))
    expect_equal(dim(decompress.res$prop),c(100,4))
    expect_equal(dim(decompress.res$sig)[[2]],4)



})
