<img align="top" style="margin-left: 1px; margin-bottom: 1px; margin-right: 1px; margin-top: 10px" src="inst/images/eCV_logo.png" width="150" height="150"/>

# Enhanced Coefficient of Variation and IDR Extensions for Reproducibility Assessment


This package provides extensions and alternative methods to measure the
reproducibility  of  omic  data  with an arbitrary number of replicates. 
It introduces an enhanced Coefficient of Variation (**eCV**)  metric to 
assess the likelihood of omic features being reproducible. Additionally, 
it offers alternatives to  the  Irreproducible  Discovery Rate (**IDR**)
calculations for multi-replicate experiments.  These tools are valuable
for analyzing high-throughput  data in genomics  and other omics fields.

## [Article](https://www.biorxiv.org/content/10.1101/2023.12.18.572208v1)

## [Vignette](https://github.com/eclipsebio/eCV/blob/main/inst/eCV_vignette.pdf)

## Install developer version.

```
# Copy and paste in R session.
if (require("remotes") == FALSE) install.packages("remotes")
install_github("eclipsebio/eCV")
```

## Load package in your R session.

```
library("eCV")
```

