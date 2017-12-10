context("plotDis")

test_that("plotDis works", {
  ranked <- rankExpr(toy_expr)
  scoredf <- singscoring(ranked, upSet = toy_up, downSet = toy_dn)
  plt1 <- plotDispersion(scoredf)
  plt2 <- plotDispersion(scoredf, annot = c(1,2))
  
  testthat::expect_true(ggplot2::is.ggplot(plt1))
  testthat::expect_true(ggplot2::is.ggplot(plt2))
  
  # annotation must have same length with samples if not NULL
  testthat::expect_error(plotDispersion(scoredf, annot = c(1)))
  testthat::expect_error(plotDispersion(scoredf, annot = c(1,2,3)))
  
})
