context("plotScoreLandscape")

test_that("plotScoreLandscape works", {
  ranked <- rankGenes(toy_expr)
  scoredf <- singscoring(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
  scoredf2 <- singscoring(ranked, upSet = toy_gs_up)
  
  
  expect_true(ggplot2::is.ggplot(plotScoreLandscape(scoredf, scoredf2)))
  expect_true(ggplot2::is.ggplot(plotScoreLandscape(
    scoredf, scoredf2,
    scorenames = c('n1', 'n2')
  )))
  expect_error(ggplot2::is.ggplot(plotScoreLandscape(
    scoredf,
    scoredf2,
    scorenames = c('n1', 'n', 'n2')
  )))
  
})
test_that("input checkings for plotScoreLandscape work", {
  ranked <- rankGenes(toy_expr)
  scoredf <- singscoring(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
  scoredf2 <- singscoring(ranked, upSet = toy_gs_up)
  
  
  expect_error((plotScoreLandscape(scoredf[1,,drop = FALSE], scoredf2)))
  rownames(scoredf2) <- c("s","t")
  expect_error((plotScoreLandscape(scoredf, scoredf2)))
})
