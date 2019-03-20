# singscore 1.3.2
* allow continous and discrete annotations for `plotDispersion` and `projectScoreLandscape`. Annotations can now be part of the score data.frame and then specified as a column name
* plot themes updated
* citation updated: cite the singscore manuscript
* number of bins for the hexbin plot in `plotLandscape` is determined from the data
* fixed bug in calculation of scores. Boundary calculation was previously done with all genes in the gene-set. It should be done with genes present in both the gene-set and the data (i.e. after filtering out those not measured in the data).
* TotalDispersion now estimated as the mean of dispersions from the up- and down-regulated gene sets instead of the sum (previous estimate divided by 2)

# singscore 1.2.2
* created a website for the package
* added ORCID IDs for authors
* added sticker for package

# singscore 0.99.9
* remove internal keyword from generateNull

# singscore 0.99.8
* vectorize getPvals
* make toy_epxr_se

# singscore 0.99.5
* represent expression data in SummarizedExperiement dataset
* function input checkings
* R code re-organised

# singscore 0.99.3
* change argument name bidirectional to knownDirection, default as TRUE

# version 0.99.2
* bidirection flag added to simpleScore() function
* optimise generateNull() function
* change GSEABase from import to Depends
