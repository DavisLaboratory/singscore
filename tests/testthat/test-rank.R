context('test generic rankExpr with different signatures')

test_that(("Test matrix input"), {
  # 
  matrixData <- as.matrix(c(1,2,5,5))
  matrixMin <- rankExpr(matrixData, tiesMethod = 'min')
  matrixDefa <- rankExpr(matrixData)
  matrixMax <- rankExpr(matrixData, tiesMethod = 'max')
  
  
  
  expect_that(matrixMin, equals(as.matrix(c(1,2,3,3))))
  expect_that(matrixDefa, equals(as.matrix(c(1,2,3,3))))
  expect_that(matrixMax, equals(as.matrix(c(1,2,4,4))))

})

test_that("Test data.frame input for rankExpr",{
  #Test data.frame input
  df <- as.data.frame(c(1,2,5,5))
  colnames(df) <- 'test'
  dfmin <- as.matrix(c(1,2,3,3))
  colnames(dfmin) <- 'test'
  dfmax <- as.matrix(c(1,2,4,4))
  colnames(dfmax) <- 'test'
  
  
  dfrMin <- rankExpr(df, tiesMethod = 'min')
  dfrDefa <- rankExpr(df)
  dfrMax <- rankExpr(df, tiesMethod = 'max')
  
  expect_that(dfrMin, equals(dfmin))
  expect_that(dfrDefa, equals(dfmin))
  expect_that(dfrMax, equals(dfmax))
})

test_that("test expressSet input for rankExpr", {
  matrixData <- as.matrix(c(1,2,5,5))
  matrixMin <- rankExpr(matrixData, tiesMethod = 'min')
  matrixDefa <- rankExpr(matrixData)
  matrixMax <- rankExpr(matrixData, tiesMethod = 'max')
  #Test S4 class 
  e <- Biobase::ExpressionSet(assayData = as.matrix(c(1,2,5,5)))
  expect_that(rankExpr(e, 'min'), is_equivalent_to(matrixMin)) 
  expect_that(rankExpr(e), is_equivalent_to(matrixDefa)) 
  expect_that(rankExpr(e, 'max'), is_equivalent_to(matrixMax)) 
})
test_that("DGEList input for rankExpr", {
  
  y <- matrix(c(5,3,2,2),ncol=2)
  dge <-  edgeR::DGEList(counts = y)
  dgeMin <- matrix(c(2,1,1,1), ncol = 2)
  dgeMax <- matrix(c(2,1,2,2), ncol = 2)
  
  expect_that(rankExpr(dge, 'min'), is_equivalent_to(dgeMin))
  expect_that(rankExpr(dge), is_equivalent_to(dgeMin))
  expect_that(rankExpr(dge, 'max'), is_equivalent_to(dgeMax))
  
})
