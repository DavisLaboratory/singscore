context("generateNull")

test_that("generateNull works", {
  ranked <- rankGenes(toy_expr)
  n_up = length(GSEABase::geneIds(toy_gs_up))
  n_down = length(GSEABase::geneIds(toy_gs_dn))
  # find out what backends can be registered on your machine
  # the first one is the default backend, and it can be changed explicitly.
  BiocParallel::registered()
 
  
  # call the permutation function to generate the empirical scores for B times.
  permuteResult = generateNull(n_up = n_up, n_down = n_down, ranked, B = 10, 
                               ncores = 1, seed = 1, 
                               useBPPARAM = BiocParallel::registered()[[1]]) 
  expect_equivalent(dim(permuteResult),c(10,ncol(ranked)))

})
