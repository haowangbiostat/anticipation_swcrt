# On Anticipation Effect in Stepped Wedge Cluster Randomized Trials
This repository contains the R code to reproduce the results presented in [On Anticipation Effect in Stepped Wedge Cluster Randomized Trials](https://github.com/haowangbiostat/anticipation_swcrt).

### Quickstart

To reproduce the results, download this repo on a machine with R, run each of R script in the [`code`](code) without modification, and then the results are saved in [`figures`](figures) or [`result`](result). All the R scripts can be run standalone. To run the R scripts, you do not need to set any pathnames; everything is relative. Only standard libraries (dplyr, geomtextpath, ggplot2, gridExtra, lme4, lmerTest, MASS,  RColorBrewer, reshape2, scales, showtext, tidyverse) are required in the R script.

### Generate illustrative figures

- Run [`figure_2.R`](code/figure_2.R) to get [`Figure 2`](figures/figure_HH.pdf) in the main article.
  - Four types of true treatment effect curves with their estimated effect curves under the HH working model.

- Run [`figure_4.R`](code/figure_4.R) to get [`Figure 4`](figures/figure_ETI.pdf) in the main article.
  - Four types of true treatment effect curves with their estimated effect curves under the ETI working model.

- Run [`figure_5.R`](code/figure_5.R) to get [`Figure 5`](figures/figure_HH-ANT.pdf) in the main article.
  - Four types of true treatment effect curves with their estimated effect curves under the HH-ANT working model.
 
### Explore coefficients under model misspecification

- Run [`figure_3.R`](code/figure_5.R) to get [`Figure 3`](figures/figure_coeff_HH.pdf) in the main article.
  - Coefficients under the HH working model when the true model is ETI.
  
- Run [`figure_8.R`](code/figure_8.R) to get [`Figure 8`](figures/figure_coeff_HH-ANT.pdf) in the main article.
  - Coefficients under the HH-ANT working model when the true model is ETI.
 
- Run [`figure_9.R`](code/figure_8.R) to get [`Figure 9`](figures/figure_coeff_HH_vs_HH-ANT.pdf) in the main article.
  - Coefficients under the HH or HH-ANT working model when the true model is ETI-ANT.
 
### Compare variances under different models

- Run [`figure_6.R`](code/figure_6.R) to get [`Figure 3`](figures/figure_variance_inflation.pdf) in the main article.
  - Contour plots of variance inflation (HH-ANT vs HH or ETI-ANT vs ETI).
 
### Compare power under different models

- Run [`figure_10.R`](code/figure_10.R) to get [`Figure 10`](figures/figure_power_ratio_delta0.01.pdf) in the main article.
- Run [`figure_11.R`](code/figure_11.R) to get [`Figure 11`](figures/figure_power_ratio_delta0.04.pdf) in the main article.
- Run [`figure_12.R`](code/figure_12.R) to get [`Figure 12`](figures/figure_power_ratio_0.2.pdf) in the main article.
- Run [`figure_13.R`](code/figure_13.R) to get [`Figure 13`](figures/figure_power_ratio_0.3.pdf) in the main article.
- Run [`figure_14.R`](code/figure_14.R) to get [`Figure 14`](figures/figure_power_ratio_0.4.pdf) in the main article.
  - Contourplots of power ratio between the treatment effect estimator from the HH-ANT and HH models.
 
- Run [`figure_6.R`](code/figure_6.R) to get [`Figure 3`](figures/figure_variance_inflation.pdf) in the main article.
  - Contour plots of variance inflation (HH-ANT vs HH or ETI-ANT vs ETI)
 
### Reproduce simulation studies

### Trial planning software
