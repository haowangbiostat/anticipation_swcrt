# On Anticipation Effect in Stepped Wedge Cluster Randomized Trials
This repository contains the R code to reproduce the results presented in [On Anticipation Effect in Stepped Wedge Cluster Randomized Trials](https://github.com/haowangbiostat/anticipation_swcrt).

### Quickstart

To reproduce the results, download this repo on a machine with R, run each of R script in the [`code`](code) without modification, and then the results are saved in [`result`](result). All the R scripts can be run standalone. To run the R scripts, you do not need to set any pathnames; everything is relative. Only standard libraries (dplyr, geomtextpath, ggplot2, gridExtra, lme4, lmerTest, MASS,  RColorBrewer, reshape2, scales, showtext, tidyverse) are required in the R script.

### Generate illustrative figures

- Run [`figure_2.R`](code/figure_2.R) to get [`Figure 2`](figures/figure_HH.pdf) in the main article
  - Four types of true treatment effect curves with their estimated effect curves under the HH working model.

- Run [`figure_4.R`](code/figure_4.R) to get [`Figure 4`](figures/figure_ETI.pdf) in the main article
  - Four types of true treatment effect curves with their estimated effect curves under the HH working model.

### Generate illustrative figures

### Reproduce simulation studies

### Trial planning software
