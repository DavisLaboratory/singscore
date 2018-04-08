context("generateNull")

test_that("generateNull works", {
  ranked <- rankGenes(toy_expr)
 
  n_down = length(GSEABase::geneIds(toy_gs_dn))
  # find out what backends can be registered on your machine
  # the first one is the default backend, and it can be changed explicitly.
  BiocParallel::registered()
 
  
  # call the permutation function to generate the empirical scores for B times.
  permuteResult = generateNull(upSet = toy_gs_up, downSet = toy_gs_dn, ranked, 
                               B = 10, ncores = 1, seed = 1, 
                               useBPPARAM = BiocParallel::registered()[[1]]) 
  expect_equivalent(dim(permuteResult),c(10,ncol(ranked)))

})
test_that("Generic features for generateNull",{
  ranked <- rankGenes(toy_expr)
  geneIdUp <- GSEABase::geneIds(toy_gs_up)
  geneIdDn <- GSEABase::geneIds(toy_gs_dn)
  
  plt <- generateNull(upSet = toy_gs_up, downSet = toy_gs_dn, ranked, B = 10, 
                      ncores = 1, seed = 1, 
                      useBPPARAM = BiocParallel::registered()[[1]]) 
  plt1 <-generateNull(upSet = geneIdUp, downSet = geneIdDn, ranked, B = 10, 
                      ncores = 1, seed = 1, 
                      useBPPARAM = BiocParallel::registered()[[1]]) 
  plt2 <- generateNull(upSet = geneIdUp, downSet = geneIdDn, ranked, B = 10, 
                       ncores = 1, seed = 1, 
                       useBPPARAM = BiocParallel::registered()[[1]]) 

  expect_equivalent(dim(plt),c(10,ncol(ranked)))
  expect_equivalent(dim(plt1),c(10,ncol(ranked)))
  expect_equivalent(dim(plt2),c(10,ncol(ranked)))
  
})