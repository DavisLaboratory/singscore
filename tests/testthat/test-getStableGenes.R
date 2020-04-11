context('test generic getStableGenes function')

test_that(("Test different types"), {
  expect_length(getStableGenes(5), 5)
  expect_length(getStableGenes(5, type = 'blood'), 5)
  expect_length(getStableGenes(5, id = 'ensembl'), 5)
  expect_warning(getStableGenes(1e5))
  expect_warning(getStableGenes(1e5, type = 'blood'))
  expect_match(getStableGenes(5, id = 'ensembl'), '^ENSG')
})