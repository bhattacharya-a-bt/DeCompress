test_that("inform set works", {

    dir = system.file('extdata',package = 'DeCompress')
    gtex = readRDS(paste0(dir,'/example_data.RDS'))
    ref = gtex$true.mixed.expression
    target = gtex$observed.mixed

    inform.set.var = findInformSet(ref,
                                   method = 'variance',
                                   n_genes = 100,
                                   n.types = 4,
                                   scree = 'cumvar')

    inform.set.lin = findInformSet(ref,
                               method = 'linearity',
                               n_genes = 100,
                               n.types = 4,
                               scree = 'cumvar')

    expect_equal(class(inform.set.var), c('matrix','array'))
    expect_equal(ncol(inform.set.var), ncol(ref))

    expect_equal(class(inform.set.lin), c('matrix','array'))
    expect_equal(ncol(inform.set.lin), ncol(ref))


})
