context('test generic-rankExpr')

test_that(("test generic rankExpr signatures"), {
  matrixDt <- as.matrix(c(1,2,5,5))
  rkin <- rankExpr(matrixDt, tiesMethod = 'min')
  rkdefault <- rankExpr(matrixDt)
  rkmx <- rankExpr(matrixDt, tiesMethod = 'max')
  expect_that(rkin,equals(as.matrix(c(1,2,3,3))))
  expect_that(rkdefault,equals(as.matrix(c(1,2,3,3))))
  expect_that(rkmx,equals(as.matrix(c(1,2,4,4))))
  
})
