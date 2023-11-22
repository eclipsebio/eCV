<img align="top" style="margin-left: 1px; margin-bottom: 1px; margin-right: 1px; margin-top: 10px" src="inst/images/eCV_logo.png" width="150" height="150"/>

# Enhanced Coefficient of Variation and IDR Extensions for Reproducibility Assessment



This package provides extensions and alternative methods to measure the
reproducibility  of  omic  data  with an arbitrary number of replicates. 
It introduces an enhanced Coefficient of Variation (**eCV**)  metric to 
assess the likelihood of omic features being reproducible. Additionally, 
it offers alternatives to  the  Irreproducible  Discovery Rate (**IDR**)
calculations for multi-replicate experiments.  These tools are valuable
for analyzing high-throughput  data in genomics  and other omics fields.

## Installing and loading the package from GitHub

```
# Copy and paste in R session.
if (require("remotes") == FALSE) install.packages("remotes")
install_github("eclipsebio/eCV")
library("eCV")
```
