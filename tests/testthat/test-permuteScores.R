context("permuteScores")

test_that("permuteScores works", {
  ranked <- rankExpr(toy_expr)
  n_up = length(GSEABase::geneIds(toy_up))
  n_down = length(GSEABase::geneIds(toy_dn))
  # find out what backends can be registered on your machine
  # the first one is the default backend, and it can be changed explicitly.
  BiocParallel::registered()
 
  
  # call the permutation function to generate the empirical scores for B times.
  permuteResult = permuteScores(n_up = n_up, n_down = n_down, ranked, B =10,
  seed = 1) 
  expect_equivalent(dim(permuteResult),c(10,ncol(ranked)))

})
