context("test-scoring")

test_that("score calculation works well with vector ids", {
  df <- as.data.frame(c(1,2,5,5))
  colnames(df) <- 'test'
  dfrMin <- rankGenes(df, tiesMethod = 'min')
  rownames(dfrMin) <- c(1,2,3,4)
  scoredfUp <- simpleScore(dfrMin, upSet = c(3,4), centerScore = TRUE)
  scoredfBoth <- simpleScore(dfrMin, upSet = c(3,4), downSet = c(1))
  
  expect_that(dim(scoredfUp), equals(c(1,2)))
  expect_that(scoredfUp, is_equivalent_to(data.frame(0.25, 0)))
  expect_that(dim(scoredfBoth), equals(c(1,6)))
  expect_that(scoredfBoth, 
              is_equivalent_to(data.frame(0.75, 0, 0.25, 0, 0.5,0)))
})

test_that("score calculation works well with GeneSet S4 object", {
  df <- as.data.frame(c(1,2,5,5))
  colnames(df) <- 'test'
  dfrMin <- rankGenes(df, tiesMethod = 'min')
  rownames(dfrMin) <- c(1,2,3,4)
  scoredfUp <- simpleScore(dfrMin, 
                           upSet = GSEABase::GeneSet(as.character(c(3,4))), 
                           centerScore = FALSE)
  scoredfBoth <- simpleScore(dfrMin, 
                             upSet = GSEABase::GeneSet(as.character(c(3,4))), 
                             downSet = GSEABase::GeneSet(as.character(c(1))),
                             centerScore = TRUE)
  expect_that(dim(scoredfUp), equals(c(1,2)))
  expect_that(scoredfUp, is_equivalent_to(data.frame(0.75, 0)))
  expect_that(dim(scoredfBoth), equals(c(1,6)))
  expect_that(scoredfBoth, 
              is_equivalent_to(data.frame(0.75, 0, 0.25, 0, 0.5,0)))
})