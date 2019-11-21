Forked version modifications:
  - Contains modified [`music_prop_scaled`](R/modified_music_prop.R) function that permits users to perform celltype batch, wherein the dataset can be subsetted by celltype into a number of expressionSets. I.e. Each expressionSet contains all samples and all genes commonly expressed (non-zero expression) across all samples, with only a few cell types from the entire dataset. This works well for large datasets that have > 50,000 cells/nuclei in total.

Multi-subject Single Cell deconvolution (MuSiC)
=============================================

`MuSiC` is a deconvolution method that utilizes cross-subject scRNA-seq to estimate cell type proportions in bulk RNA-seq data.
![MuSiC\_pipeline](FigureMethod.jpg)

How to cite `MuSiC`
-------------------
Please cite the following publication:

> *Bulk tissue cell type deconvolution with multi-subject single-cell expression reference*<br />
> <small>X. Wang, J. Park, K. Susztak, N.R. Zhang, M. Li<br /></small>
> Nature Communications. 2019 Jan 22 [https://doi.org/10.1038/s41467-018-08023-x](https://doi.org/10.1038/s41467-018-08023-x) 

Installation
------------

``` r
# install devtools if necessary
install.packages('devtools')

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# load
library(MuSiC)
```

More Information
-----------------
Please see [Tutorial](http://xuranw.github.io/MuSiC/articles/MuSiC.html).
