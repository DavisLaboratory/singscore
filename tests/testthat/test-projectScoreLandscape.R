context("projectScoreLandscape")

test_that("projectScoreLandscape works", {
  ranked <- rankGenes(toy_expr)
  scoredf1 <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
  scoredf2 <- simpleScore(ranked, upSet = toy_gs_up)
  psl <- plotScoreLandscape(scoredf1, scoredf2)
  
  expect_true(ggplot2::is.ggplot(
    suppressWarnings(projectScoreLandscape(psl,scoredf1, scoredf2))))
  expect_true(ggplot2::is.ggplot(
    suppressWarnings(projectScoreLandscape(psl,scoredf1, scoredf2, c(1)))))
  
  #sampleLabels mush have same number of elements with samples
  expect_error(ggplot2::is.ggplot(
    suppressWarnings(projectScoreLandscape(psl,scoredf1, 
                                           scoredf2, c(1:5),
                                           sampleLabels = c('l1','l2','l3',
                                                            'l4','l5')))))
  expect_error(ggplot2::is.ggplot(
    suppressWarnings(projectScoreLandscape(psl,scoredf1, 
                                           scoredf2,
                                           sampleLabels = c('l1','l2',
                                                            'l4','l5')))))
  expect_error(ggplot2::is.ggplot(
    suppressWarnings(projectScoreLandscape(psl,scoredf1, 
                                           scoredf2,
                                           annot = c('a1','a2','a2')))))
  
  
})
