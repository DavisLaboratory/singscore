context("testing plotRankDensity function ")

test_that("plotRankDensity is plotting", {
  ranked <- rankExpr(toy_expr)
  plt <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_up)
  plt2 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_up,
                         downSet = toy_dn)
  expect_true(ggplot2::is.ggplot(plt))
  expect_true(ggplot2::is.ggplot(plt2))
 # expect_error(plotRankDensity(ranked[,2], upSet = toy_up))
})

test_that("Generic features for plotRankDensity",{
  ranked <- rankExpr(toy_expr)
  geneIdUp <- GSEABase::geneIds(toy_up)
  geneIdDn <- GSEABase::geneIds(toy_dn)
  
  plt <- plotRankDensity(ranked[,2,drop = FALSE], upSet = geneIdUp)
  plt1 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_up)
  plt2 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = geneIdUp,
                          downSet = toy_dn)
  plt3 <- plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_up,
                          downSet = geneIdDn)
  expect_true(ggplot2::is.ggplot(plt))
  expect_true(ggplot2::is.ggplot(plt2))
  expect_true(ggplot2::is.ggplot(plt3))
  expect_that(plt2, is_equivalent_to(plt3))
  expect_that(plt, is_equivalent_to(plt1))
  
})