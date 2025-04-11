
<!-- README.md is generated from README.Rmd. Please edit that file -->

# devgru <img src="man/figures/devgruLogo.png" align="right" width="250" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/pblaney/devgru)](https://github.com/pblaney/devgru/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/pblaney/devgru)](https://github.com/pblaney/devgru/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

A developmental `R` environment for a suite of `GenomicRanges` utilities

### Installation

You can install the current development version of `devgru` with:

``` r
devtools::install_github("pblaney/devgru")
#BiocManager::install("pblaney/devgru")
```

### Suite of Packages

The current suite of utilities for `devgru` can be quickly installed, if
needed, and loaded with:

``` r
kit_loadout()
```

#### GenomicRanges Core

- [`BSgenome.Hsapiens.UCSC.hg38`](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html)
- [`GenomicRanges`](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- [`GenomeInfoDb`](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html)
- [`data.table`](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
- [`mskilab-org/gUtils`](https://github.com/mskilab-org/gUtils)
- [`VariantAnnotation`](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
- [`rtracklayer`](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [`Biostrings`](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [`S4Vectors`](https://bioconductor.org/packages/release/bioc/html/S4Vectors.html)

#### Utility Core

- [`dplyr`](https://dplyr.tidyverse.org)
- [`stringr`](https://stringr.tidyverse.org)
- [`readr`](https://readr.tidyverse.org)
- [`ggplot2`](https://ggplot2.tidyverse.org)
- [`ggsci`](https://nanx.me/ggsci/)
- [`paletter`](https://emilhvitfeldt.github.io/paletteer/)
- [`scico`](https://github.com/thomasp85/scico)
- [`flextable`](https://ardata-fr.github.io/flextable-book/index.html)
- [`mclust`](https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html)
- [`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
- [`doParallel`](https://cran.r-project.org/web/packages/doParallel/doParallel.pdf)
- [`foreach`](https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html)

### Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.14/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://pblaney.github.io/devgru) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.14/biocthis)*.
