context("plotDis")

test_that("plotDispersion checks works", {
  ranked <- rankGenes(toy_expr_se)
  scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
  plt1 <- plotDispersion(scoredf)
  plt2 <- plotDispersion(scoredf, annot = c(1,2))

  testthat::expect_true(ggplot2::is_ggplot(plt1))
  testthat::expect_true(ggplot2::is_ggplot(plt2))

  # annotation must have same length with samples if not NULL
  testthat::expect_error(plotDispersion(scoredf, annot = c(1)))
  testthat::expect_error(plotDispersion(scoredf, annot = c(1,2,3)))
  testthat::expect_error(plotDispersion(scoredf, textSize = "ss"))
})

test_that('plotDispersion colours work for single sigs', {
	#score
	scoredf = scoredf_ccle_epi

	#annotations
	discrete_annot = as.factor(sample.int(3, nrow(scoredf), replace = TRUE))
	continous_annot = rnorm(nrow(scoredf))
	char_annot = LETTERS[discrete_annot]

	p1 = plotDispersion(scoredf)
	testthat::expect_true(ggplot2::is_ggplot(p1))

	p2 = plotDispersion(scoredf, annot = discrete_annot)
	testthat::expect_true(ggplot2::is_ggplot(p2))

	p3 = plotDispersion(scoredf, annot = continous_annot)
	testthat::expect_true(ggplot2::is_ggplot(p3))

	p4 = plotDispersion(scoredf, annot = char_annot)
	testthat::expect_true(ggplot2::is_ggplot(p4))

	#column annotation
	scoredf$MyAnnot = char_annot
	p5 = plotDispersion(scoredf, annot = 'MyAnnot')
	testthat::expect_true(ggplot2::is_ggplot(p5))
})

test_that('plotDispersion colours work for up/dn sigs', {
	#score
	scoredf = as.data.frame(matrix(runif(60), 10))
	cnames = paste0(rep(c('Total', 'Up', 'Down'), each = 2), c('Score', 'Dispersion'))
	colnames(scoredf) = cnames

	#annotations
	discrete_annot = as.factor(sample.int(3, nrow(scoredf), replace = TRUE))
	continous_annot = rnorm(nrow(scoredf))
	char_annot = LETTERS[discrete_annot]

	p1 = plotDispersion(scoredf)
	testthat::expect_true(ggplot2::is_ggplot(p1))

	p2 = plotDispersion(scoredf, annot = discrete_annot)
	testthat::expect_true(ggplot2::is_ggplot(p2))

	p3 = plotDispersion(scoredf, annot = continous_annot)
	testthat::expect_true(ggplot2::is_ggplot(p3))

	p4 = plotDispersion(scoredf, annot = char_annot)
	testthat::expect_true(ggplot2::is_ggplot(p4))

	#column annotation
	scoredf$MyAnnot = char_annot
	p5 = plotDispersion(scoredf, annot = 'MyAnnot')
	testthat::expect_true(ggplot2::is_ggplot(p5))
})

test_that('plotDispersion interactive for single sigs', {
	#score
	scoredf = scoredf_ccle_epi

	#annotations
	discrete_annot = as.factor(sample.int(3, nrow(scoredf), replace = TRUE))
	continous_annot = rnorm(nrow(scoredf))
	char_annot = LETTERS[discrete_annot]

	p1 = plotDispersion(scoredf, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p1))

	p2 = plotDispersion(scoredf, annot = discrete_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p2))

	p3 = plotDispersion(scoredf, annot = continous_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p3))

	p4 = plotDispersion(scoredf, annot = char_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p4))

	#column annotation
	scoredf$MyAnnot = char_annot
	p5 = plotDispersion(scoredf, annot = 'MyAnnot', isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p5))
})

test_that('plotDispersion interactive for up/dn sigs', {
	#score
	scoredf = as.data.frame(matrix(runif(60), 10))
	cnames = paste0(rep(c('Total', 'Up', 'Down'), each = 2), c('Score', 'Dispersion'))
	colnames(scoredf) = cnames

	#annotations
	discrete_annot = as.factor(sample.int(3, nrow(scoredf), replace = TRUE))
	continous_annot = rnorm(nrow(scoredf))
	char_annot = LETTERS[discrete_annot]

	p1 = plotDispersion(scoredf, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p1))

	p2 = plotDispersion(scoredf, annot = discrete_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p2))

	p3 = plotDispersion(scoredf, annot = continous_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p3))

	p4 = plotDispersion(scoredf, annot = char_annot, isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p4))

	#column annotation
	scoredf$MyAnnot = char_annot
	p5 = plotDispersion(scoredf, annot = 'MyAnnot', isInteractive = TRUE)
	testthat::expect_true('plotly' %in% class(p5))
})

