
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ToDcallR

<!-- badges: start -->

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Project](https://img.shields.io/badge/Master--PhD-Thesis-blue)
![Version](https://img.shields.io/badge/Version-0.1.0-red) [![License:
MIT](https://cdn.prod.website-files.com/5e0f1144930a8bc8aace526c/65dd9eb5aaca434fac4f1c34_License-MIT-blue.svg)](/LICENSE)
<!-- badges: end -->

ToDcallR is an R package developed to call Translation on Demand (ToD)
candidate genes in coupled temporal bulk RNA-Seq and proteomics
datasets. The functions of this package correspond to the developed
workflow of one of Luke Brandwood’s PhD projects and my (Elias Schwall)
Master Thesis project (Translation on Demand as a post-transcriptional
regulation mechanism of embryonic stem cell differentiation) at the
[Beyer
Lab](https://www.cecad.uni-koeln.de/research/principal-investigators/full-members/andreas-beyer),
[CECAD](https://www.cecad.uni-koeln.de/home), [University of
Cologne](https://www.uni-koeln.de/en/).

## Translation on Demand

We define ToD as an increase in protein abundance due to enhanced
translation, occurring without any changes in mRNA levels of the
corresponding gene. In other words, ToD is a regulatory mechanism where
translation efficiency is selectively increased for transcripts of
specific genes. This upregulation allows for rapid adjustments in
protein levels, supporting cellular responses that need immediate
protein synthesis without requiring new mRNA transcription. Such
translation-driven increases could be particularly useful in situations
where quick adaptations are necessary, as it enables cells to respond to
changes in the environment or developmental cues by rapidly elevating
the protein output from existing mRNA pools.

## Installation

You can install the development version of ToDcallR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("eliasschwall/ToDcallR")
```

## How to use ToDcallR

Before using ToDcallR you should make sure that your datasets fit the
expected format:

- The transcriptomic and proteomic dataset need to have or be restricted
  to the same number of time points
- The `column names` of the dataframes should be corresponding to the
  time point in hours e.g.: `0h`,`0.5h`,`1h` etc.
- If one or both datasets have replicates indicate them by `_rep1` etc.
  e.g.: `0h_rep1`,`0h_rep2`,`0.5h_rep1`,`0.5h_rep2`

> [!IMPORTANT]
> Some thoughts on normalization: Although ultimately the form of normalization is up to you and your specific data modality, on temporal RNA-Seq or proteomics data, it is recommended to use a method that leverages information from adjacent samples, like gaussian process regression (GPR) or LOESS. This is even more important if you have a low number of replicates. It is up to you to transform the data into a form that lets you derive meaningful conclusions. All ToDCallR does is compare log fold change (LFC) values to infer how protein abundance changes in relation to transcript abundance changes for a particular gene.

``` r
library(ToDcallR)
## basic example code
```
