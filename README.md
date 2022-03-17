
singscore <img src="man/figures/logo.png" align="right"  height="140" width="120" alt="logo"/>
==============================================================================================

[![R-CMD-check](https://github.com/DavisLaboratory/singscore/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/DavisLaboratory/singscore/actions)
[![codecov](https://codecov.io/gh/DavisLaboratory/vissE/branch/main/graph/badge.svg?token=8JHZB1GN26)](https://codecov.io/gh/DavisLaboratory/vissE)
[![BioC status](https://bioconductor.org/shields/years-in-bioc/singscore.svg)](https://bioconductor.org/packages/singscore/)

Overview
--------

'singscore' is an R/Bioconductor package which implements the simple single-sample gene-set (or gene-signature) scoring method proposed by Foroutan *et al.* (2018). It uses rank-based statistics to analyze each sample's gene expression profile and scores the expression activities of gene sets at a single-sample level.

Getting Started
---------------

These instructions will get you to install the package up and running on your local machine. If you experience any issues, please raise a GitHub issue at <https://github.com/DavisLaboratory/singscore/issues>.

    # build_vignettes = TRUE to build vignettes upon installation
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("singscore")

Documentation
-------------

The package comes with a vignette showing how the different functions in the package can be used to perform a gene-set enrichment analysis on a single sample level. Pre-built vignettes can be accessed via [Bioconductor](https://bioconductor.org/packages/release/bioc/vignettes/singscore/inst/doc/singscore.html) or [the GitHub IO page](https://davislaboratory.github.io/singscore/articles/singscore.html).

A detailed workflow on using singscore to analyse the TCGA acute myeloid leukemia (AML) data is available as part of the [SingscoreAMLMutations](https://doi.org/10.18129/B9.bioc.SingscoreAMLMutations) R/Bioconductor package.

References
----------

Foroutan, Momeneh, Dharmesh D Bhuva, Ruqian Lyu, Kristy Horan, Joseph Cursons, and Melissa J Davis. 2018. “Single Sample Scoring of Molecular Phenotypes.” BMC Bioinformatics 19 (1). BioMed Central: 404. doi: [10.1186/s12859-018-2435-4](https://doi.org/10.1186/s12859-018-2435-4).
