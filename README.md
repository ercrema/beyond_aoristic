# Data and R scripts for the paper 'A Bayesian alternative for Aoristic analyses in Archaeology'

This repository contains data and scripts used in the following paper:

Crema, E.R. (2023) A Bayesian alternative for Aoristic analyses in Archaeology.

The repository is organised into four directories: _figures_, _results_, _scripts_, and _src_. 
The _scripts_ directory contains R scripts for each of the four simulation experiments, _src_ contains utility functions, _results_ contains R image files with the analyses output, and _figures_ contains all figures in the manuscript and the R script required to generate them. 

## File Structure

#### figures
* `figure1.png` ~ `figure6.png` ... PNG version of manuscript figures. Generated using the command `Rscript figures.R png`. TIFF version can be generated with `Rscript figures.R tiff`.
* `figures.R` ... generates all figures for the manuscript.
  
#### results
* `exo01_res.RData` ... R image file containing results of experiment 1, generated running `experiment01.R`.
* `exo02a_res.RData` ... R image file containing results of experiment 2 (100 yrs resolution), generated running `experiment02a.R`.
* `exo02b_res.RData` ... R image file containing results of experiment 2 (10 yrs resolution), generated running `experiment02b.R`.
* `exo03_res.RData` ... R image file containing results of experiment 3, generated running `experiment03.R`.
* `exo04_res.RData` ... R image file containing results of experiment 4, generated running `experiment04.R`.
* `figure1_res.RData` ... R image containing objects required for figure 1.
  
#### scripts
* `experiment01.R` ... R script for running experiment #1
* `experiment02a.R` ... R script for running experiment #2 (100yrs resolution).
* `experiment02b.R` ... R script for running experiment #2 (10yrs resolution).
* `experiment03.R` ... R script for running experiment #3
* `experiment04.R` ... R script for running experiment #4
  

#### src
* `diristick.R` ... Broken stick algorithm for generating random archaeological periodisation using the Dirichlet distribution. 
* `mcsim.R` ... Monte-Carlo simulation for aoristic analyses.
* `phaserect.R` ... displays periodisations interval on existing plots.
* `randtimespan.R` ... assigns random time-spans to events.
* `ribbon.R` ... displays interval-based ribbons on existing plots.
* `time2phase.R` ... assign periods/phases to calendar dates based on user-defined intervals.
* `unifdisc.R` ... discrete uniform distribution functions.


## R Session Info
```
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] nimbleCarbon_0.2.4 emdbook_1.3.13     dplyr_1.1.2        baorista_0.0.6    
[5] nimble_1.0.1       here_1.0.1        

loaded via a namespace (and not attached):
 [1] xfun_0.40              spatstat.sparse_3.0-2  lattice_0.22-5        
 [4] numDeriv_2016.8-1.1    vctrs_0.6.2            tools_4.3.1           
 [7] doSNOW_1.0.20          spatstat.utils_3.0-3   generics_0.1.3        
[10] goftest_1.2-3          stats4_4.3.1           parallel_4.3.1        
[13] tibble_3.2.1           proxy_0.4-27           fansi_1.0.4           
[16] pkgconfig_2.0.3        spatstat_3.0-6         Matrix_1.6-1.1        
[19] KernSmooth_2.23-22     lifecycle_1.0.3        compiler_4.3.1        
[22] deldir_1.0-9           spatstat.linnet_3.1-1  codetools_0.2-19      
[25] spatstat.explore_3.2-1 snow_0.4-4             class_7.3-22          
[28] pracma_2.4.2           pillar_1.9.0           MASS_7.3-60           
[31] classInt_0.4-9         spatstat.model_3.2-4   iterators_1.0.14      
[34] rpart_4.1.21           abind_1.4-5            foreach_1.5.2         
[37] nlme_3.1-163           spatstat.geom_3.2-4    tidyselect_1.2.0      
[40] bdsmatrix_1.3-6        mvtnorm_1.2-2          sf_1.0-14             
[43] splines_4.3.1          polyclip_1.10-4        rprojroot_2.0.3       
[46] grid_4.3.1             cli_3.6.1              magrittr_2.0.3        
[49] utf8_1.2.3             e1071_1.7-13           spatstat.data_3.0-1   
[52] tensor_1.5             igraph_1.5.1           coda_0.19-4           
[55] knitr_1.43             rcarbon_1.5.1          bbmle_1.0.25          
[58] mgcv_1.9-0             rlang_1.1.1            Rcpp_1.0.11           
[61] spatstat.random_3.1-5  glue_1.6.2             DBI_1.1.3             
[64] R6_2.5.1               plyr_1.8.8             units_0.8-3      
```
Please note that the R package _baorista_ is not available on CRAN yet, but the latest version can be installed directly from GitHub with the following command:
```
library(devtools)
install_github('ercrema/baorista')
```
## Funding
 * ERC Starting Grant Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER) (Project N. 801953, PI: E. Crema).
 * Philip Leverhulme Prize (#PLP-2019–304 Awarded to: E.Crema)

## Licence
CC-BY 3.0
