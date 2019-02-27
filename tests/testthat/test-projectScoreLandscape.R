context("projectScoreLandscape")

test_that("projectScoreLandscape works", {
  ranked <- rankGenes(toy_expr_se)
  scoredf1 <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
  scoredf2 <- simpleScore(ranked, upSet = toy_gs_up)
  psl <- plotScoreLandscape(scoredf1, scoredf2)

  expect_true(ggplot2::is.ggplot(
    suppressWarnings(projectScoreLandscape(psl,scoredf1, scoredf2))))
  expect_error(ggplot2::is.ggplot(
    suppressWarnings(projectScoreLandscape(psl,scoredf1, scoredf2, c(1)))))

  #sampleLabels must have same number of elements with samples
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

test_that('projectScoreLandscape colours work for different annotations', {
	#score
	p0 = plotScoreLandscape(scoredf_tcga_epi, scoredf_tcga_mes)

	#annotations
	discrete_annot = as.factor(sample.int(3, nrow(scoredf_ccle_epi), replace = TRUE))
	continous_annot = rnorm(nrow(scoredf_ccle_epi))
	char_annot = LETTERS[discrete_annot]

	p1 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes)
	testthat::expect_true(ggplot2::is.ggplot(p1))

	p2 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes, discrete_annot)
	testthat::expect_true(ggplot2::is.ggplot(p2))

	p3 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes, continous_annot)
	testthat::expect_true(ggplot2::is.ggplot(p3))

	p4 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes, char_annot)
	testthat::expect_true(ggplot2::is.ggplot(p4))

	#column annotation
	scoredf = scoredf_ccle_epi
	scoredf$MyAnnot = char_annot
	p5 = projectScoreLandscape(p0, scoredf, scoredf_ccle_mes, 'MyAnnot')
	testthat::expect_true(ggplot2::is.ggplot(p5))
})

test_that('projectScoreLandscape works for interactive plots', {
	#score
	p0 = plotScoreLandscape(scoredf_tcga_epi, scoredf_tcga_mes)

	#annotations
	discrete_annot = as.factor(sample.int(3, nrow(scoredf_ccle_epi), replace = TRUE))
	continous_annot = rnorm(nrow(scoredf_ccle_epi))
	char_annot = LETTERS[discrete_annot]

	p1 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p1))

	p2 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes, discrete_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p2))

	p3 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes, continous_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p3))

	p4 = projectScoreLandscape(p0, scoredf_ccle_epi, scoredf_ccle_mes, char_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p4))

	#column annotation
	scoredf = scoredf_ccle_epi
	scoredf$MyAnnot = char_annot
	p5 = projectScoreLandscape(p0, scoredf, scoredf_ccle_mes, 'MyAnnot', isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p5))
})