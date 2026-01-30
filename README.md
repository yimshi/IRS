# IRS: Iterative Reference Selection for Microbiome Data

The `IRS` package provides a robust statistical framework for identifying stable reference taxa in microbiome count data. In compositionality-affected data, selecting an appropriate reference set is crucial for unbiased normalization and valid downstream differential abundance testing.



## Methodology

IRS employs a dual-filter iterative approach to ensure the selected reference taxa are independent of the study's primary predictor:

* **Non-parametric Filter**: Uses Kendall's Tau to assess monotonic relationships on normalized abundances, ensuring robustness against outliers and heavy tails.
* **Robust Parametric Filter**: Implements a Poisson Generalized Linear Model (GLM) with **sandwich (HC3) variance estimators** to handle overdispersion and model misspecification common in sequencing counts.

## Installation

You can install the development version of IRS from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("yimshi/IRS")
