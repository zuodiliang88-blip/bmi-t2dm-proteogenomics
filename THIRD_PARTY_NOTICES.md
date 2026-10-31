# Third-party code

`CARMA_strict.R` is the study authors' pre-existing modification of CARMA, originally authored by Zikun Yang, Chen Wang and Iuliana Ionita-Laza. The installed CARMA 1.0 DESCRIPTION records `GPL (>= 2)`; its accompanying LICENSE contains GPL version 3. This distribution includes that licence verbatim as `LICENSE-CARMA` and distributes the modified CARMA source under GPL-3.0-or-later.

The study modification enforces the maximum number of causal variants during model-space exploration. It predates this repository; packaging on 2026-10-31 only removed comments/whitespace and added attribution, with executable expressions unchanged. Callers explicitly request a maximum of five causal variants and enable outlier detection.

Reference: Yang Z, Wang C, Liu L, et al. CARMA is a new Bayesian model for fine-mapping in genome-wide association meta-analyses. *Nature Genetics* 55, 1057–1065 (2023). https://doi.org/10.1038/s41588-023-01392-0

`moloc_typed.R` copies and modifies a function from the installed moloc package at runtime to preserve the explicit trait types supplied by the study. It does not bundle the package implementation. This pre-existing study patch is retained without changing its executable expressions.
