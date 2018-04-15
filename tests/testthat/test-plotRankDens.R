context("testing plotRankDensity function ")

test_that("plotRankDensity is plotting", {
  ranked <- rankGenes(toy_expr_se)
  plt <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_gs_up)
  plt2 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_gs_up,
                         downSet = toy_gs_dn)
  expect_true(ggplot2::is.ggplot(plt))
  expect_true(ggplot2::is.ggplot(plt2))
 # expect_error(plotRankDensity(ranked[,2], upSet = toy_gs_up))
})

test_that("Generic features for plotRankDensity",{
  ranked <- rankGenes(toy_expr_se)
  geneIdUp <- GSEABase::geneIds(toy_gs_up)
  geneIdDn <- GSEABase::geneIds(toy_gs_dn)
  
  plt <- plotRankDensity(ranked[,2,drop = FALSE], upSet = geneIdUp)
  plt1 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_gs_up)
  plt2 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = geneIdUp,
                          downSet = geneIdDn)
  plt3 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_gs_up,
                          downSet = toy_gs_dn)
  expect_true(ggplot2::is.ggplot(plt))
  expect_true(ggplot2::is.ggplot(plt2))
  expect_true(ggplot2::is.ggplot(plt3))
  expect_that(plt2, is_equivalent_to(plt3))
  expect_that(plt, is_equivalent_to(plt1))
  
})

test_that("input checkings work for plotRankDensity",{
  ranked <- rankGenes(toy_expr_se)
  geneIdUp <- GSEABase::geneIds(toy_gs_up)
  geneIdDn <- GSEABase::geneIds(toy_gs_dn)

  expect_error(plotRankDensity(ranked[,2,drop = FALSE], upSet = geneIdUp, 
                               textSize = "ss"))
  expect_error(plotRankDensity(ranked[,2,drop = FALSE], upSet = geneIdUp, 
                               isInteractive = "KK"))
  
})