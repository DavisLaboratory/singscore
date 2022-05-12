
# singscore <img src="man/figures/logo.png" align="right"  height="140" width="120" alt="logo"/>

## Overview

‘singscore’ is an R/Bioconductor package which implements the simple
single-sample gene-set (or gene-signature) scoring method proposed by
Foroutan *et al.* (2018) and Bhuva *et al.* (2020). It uses rank-based
statistics to analyze each sample’s gene expression profile and scores
the expression activities of gene sets at a single-sample level.

Additional packages we have developed can enhance the singscore
workflow:

1.  [`msigdb`](https://www.bioconductor.org/packages/release/data/experiment/html/msigdb.html) -
    A package that provides gene-sets from the molecular signatures
    database (MSigDB) as a `GeneSetCollection` object that is compatible
    with `singscore`.
2.  [`vissE`](https://www.bioconductor.org/packages/release/bioc/html/vissE.html) -
    A package that can summarise and aid in the interpretation of a list
    of significant gene-sets identified by `singscore` (see
    [tutorial](https://davislaboratory.github.io/GenesetAnalysisWorkflow/)).
3.  [`emtdata`](https://www.bioconductor.org/packages/release/data/experiment/html/emtdata.html) -
    The full EMT dataset used in this tutorial (with additional EMT
    related datasets).

We have also published and made openly available the extensive tutorials
below that demonstrate the variety of ways in which `singscore` can be
used to gain a better functional understanding of molecular data:

1.  [Using singscore to predict mutation status in acute myeloid
    leukemia from transcriptomic
    signatures](https://f1000research.com/articles/8-776).
2.  [Gene-set enrichment analysis
    workshop](https://davislaboratory.github.io/GenesetAnalysisWorkflow/) -
    available through the
    [Orchestra](http://app.orchestra.cancerdatasci.org/) platform
    (search “WEHI Masterclass Day 4: Functional Analysis, single sample
    gene set analysis”).

## Getting Started

These instructions will get you to install the package up and running on
your local machine. If you experience any issues, please raise a GitHub
issue at <https://github.com/DavisLaboratory/singscore/issues>.

    # build_vignettes = TRUE to build vignettes upon installation
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("singscore", version = "3.8")

## Documentation

The package comes with a vignette showing how the different functions in
the package can be used to perform a gene-set enrichment analysis on a
single sample level. Pre-built vignettes can be accessed via
[Bioconductor](https://bioconductor.org/packages/release/bioc/vignettes/singscore/inst/doc/singscore.html)
or [the GitHub IO
page](https://davislaboratory.github.io/singscore/articles/singscore.html).

## References

Foroutan M, Bhuva D, Lyu R, Horan K, Cursons J, Davis M (2018). “Single
sample scoring of molecular phenotypes.” *BMC bioinformatics*, *19*(1),
404. doi:
[10.1186/s12859-018-2435-4](https://doi.org/10.1186/s12859-018-2435-4).

Bhuva D, Cursons J, Davis M (2020). “Stable gene expression for
normalisation and single-sample scoring.” *Nucleic Acids Research*,
*48*(19), e113. doi:
[10.1093/nar/gkaa802](https://doi.org/10.1093/nar/gkaa802).
