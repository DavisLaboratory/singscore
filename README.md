# singscore <img src="man/figures/logo.png" align="right"  height="140" width="120" alt="logo"/>

## Overview

'singscore' is a gene-set scoring R package which implements a simple single-sample gene-set (or gene-signature) scoring method. It uses rank-based statistics to analyze the each sample's gene expression profile and scores the expression activities of gene sets at a single-sample level.


## Getting Started

These instructions will get you to install the package up and running on your local machine. Please be aware that the package is still under development and we aim to submit the release version to [Bioconductor](https://www.bioconductor.org). If you experience any issues, please let us know via email 
<ul>
 <li> Ruqian Lyu lyu.r@wehi.edu.au </li>
 <li> Momeneh Foroutan foroutan.m@wehi.edu.au</li> 
 <li> Dharmesh Bhuva bhuva.d@wehi.edu.au </li>
 <li> Joseph Cursons cursons.j@wehi.edu.au</li>
 <li> Melissa Davis davis.m@wehi.edu.au</li>
</ul>

### Prerequisites
To install the package from git hub, install the package 'devtools' first and then use the function `install_github` to install 'singscore' by running the following script


### Installing

```
install.packages('devtools')
# build_vignettes = TRUE to build vignettes upon installation
devtools::install_github('DavisLaboratory/singscore', build_vignettes = TRUE)
library(singscore)
```


## Running singscoring

```
# The example data sets
# We have included several example datasets in the package to illustrate the 
# usage of scoring and visulisation functions

# an SummarizedExperiment object containing expression dataset 
str(tgfb_expr_10_se)
# Get the counts by 
SummarizedExperiment::assay(tgfb_expr_10_se)
# up gene set
str(tgfb_gs_up)
# down gene set
str(tgfb_gs_dn)

# rank the expression matrix first
rankedData <- rankGenes(tgfb_expr_10_se)

# Call simpleScore function to score each individual sample
scoredf <- simpleScore(rankedData, tgfb_gs_up,tgfb_gs_dn, centerScore = TRUE)

scoredf
##               TotalScore TotalDispersion    UpScore UpDispersion
## D_Ctrl_R1   -0.088097993        5734.697 0.06096415     3119.390
## D_TGFb_R1    0.286994210        4435.939 0.24931565     2352.886
## D_Ctrl_R2   -0.098964086        5722.836 0.06841242     3129.769
## D_TGFb_R2    0.270721958        4757.663 0.25035661     2470.012
## Hes_Ctrl_R1 -0.002084788        5492.292 0.08046490     3134.216
## Hes_TGFb_R1  0.176122839        5195.030 0.22894035     2416.638
## Hes_Ctrl_R2  0.016883867        5401.112 0.08817828     3138.664
## Hes_TGFb_R2  0.188466953        4910.371 0.23895473     2324.717
## Hil_Ctrl_R1 -0.061991164        6078.660 0.08314254     3553.792
## Hil_Ctrl_R2 -0.064937366        5918.539 0.07433863     3396.637
##               DownScore DownDispersion
## D_Ctrl_R1   -0.14906214       2615.306
## D_TGFb_R1    0.03767856       2083.053
## D_Ctrl_R2   -0.16737650       2593.067
## D_TGFb_R2    0.02036534       2287.652
## Hes_Ctrl_R1 -0.08254969       2358.075
## Hes_TGFb_R1 -0.05281751       2778.392
## Hes_Ctrl_R2 -0.07129441       2262.448
## Hes_TGFb_R2 -0.05048778       2585.654
## Hil_Ctrl_R1 -0.14513371       2524.868
## Hil_Ctrl_R2 -0.13927600       2521.903

```



## Produced rank density plot
`plotRankDensity()` takes a single column data frame, which is a subset of the ranked data obtained from `rankGenes()`function, and gene sets, and it returns plots visualising the density and the rugs of the ranks.
```
plotRankDensity(rankedData[,2,drop = FALSE], upSet = tgfb_gs_up, 
                downSet = tgfb_gs_dn, isInteractive = FALSE)

```
Which can produce 

![](https://user-images.githubusercontent.com/12887308/37870731-5583c39e-3028-11e8-9ddb-5d197c05d55a.png)

## Plot Score Landscape

plotScoreLandscape plots the scores of the samples against two different gene signatures in a landscape for exploring the replationships of gene signatures.

```
data("scoredf_tcga_epi")
data("scoredf_tcga_mes")
plotScoreLandscape(scoredf_tcga_epi, scoredf_tcga_mes, 
                   scorenames = c('tcga-EPI','tcga-MES'), isInteractive = FALSE)

```
![](https://user-images.githubusercontent.com/12887308/37870753-d12692ba-3028-11e8-9aaf-37bceef5722b.png)

**For more examples, please refer to the vigntte included in the package by**
`browseVignettes('singscore')`

### Linked Article
[Defining a landscape of molecular phenotypes using a simple single sample scoring method](https://www.biorxiv.org/content/early/2017/12/08/231217 )
