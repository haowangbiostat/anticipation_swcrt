# On Anticipation Effect in Stepped Wedge Cluster Randomized Trials
This repository contains the R code to reproduce the results presented in [On Anticipation Effect in Stepped Wedge Cluster Randomized Trials](https://github.com/haowangbiostat/anticipation_swcrt).

### Quickstart

To reproduce the results, download this repo on a machine with R, run each of R script in the [`code`](code) without modification, and then the results are saved in [`figures`](figures) or [`result`](result). All the R scripts can be run standalone. To run the R scripts, you do not need to set any pathnames; everything is relative. Only standard libraries (dplyr, geomtextpath, ggplot2, gridExtra, lme4, lmerTest, MASS,  RColorBrewer, reshape2, scales, showtext, tidyverse) are required in the R script.

### Generate Illustrative Figures

- Run [`figure_2.R`](code/figure_2.R) to get [`Figure 2`](figures/figure_HH.pdf) in the main article.
  - Four types of true treatment effect curves with their estimated effect curves under the HH working model.
  - True model is HH-ANT or ETI-ANT.

- Run [`figure_4.R`](code/figure_4.R) to get [`Figure 4`](figures/figure_ETI.pdf) in the main article.
  - Four types of true treatment effect curves with their estimated effect curves under the ETI working model.
  - True model is HH-ANT or ETI-ANT.

- Run [`figure_5.R`](code/figure_5.R) to get [`Figure 5`](figures/figure_HH-ANT.pdf) in the main article.
  - Four types of true treatment effect curves with their estimated effect curves under the HH or HH-ANT working model.
  - True model is HH, HH-ANT, ETI, or ETI-ANT.
 
- Run [`figure_18.R`](code/figure_18.R) to get [`Figure 18`](figures/figure_HH_higher_order.pdf) in the main article.
  - Constant (true) treatment effect curves with their estimated effect curves under the HH working model.
  - True model is HH-ANT with a higher-order anticipation effect.
 
### Explore Coefficients Under Model Misspecification

- Run [`figure_3.R`](code/figure_3.R) to get [`Figure 3`](figures/figure_coeff_HH.pdf) in the main article.
  - Coefficients under the HH working model when the true model is ETI.
  
- Run [`figure_8.R`](code/figure_8.R) to get [`Figure 8`](figures/figure_coeff_HH-ANT.pdf) in the main article.
  - Coefficients under the HH-ANT working model when the true model is ETI.
 
- Run [`figure_9.R`](code/figure_9.R) to get [`Figure 9`](figures/figure_coeff_HH_vs_HH-ANT.pdf) in the main article.
  - Coefficients under the HH or HH-ANT working model when the true model is ETI-ANT.
 
- Run [`figure_17.R`](code/figure_17.R) to get [`Figure 17`](figures/figure_coeff_HH_higher_order.pdf) in the main article.
  - Coefficients under the HH working model when the true model is HH-ANT with a higher-order anticipation effect.
 
### Compare Variances Under Different Models

- Run [`figure_6.R`](code/figure_6.R) to get [`Figure 6`](figures/figure_variance_inflation.pdf) in the main article.
  - Contour plots of variance inflation (HH-ANT vs HH or ETI-ANT vs ETI).
 
### Compare Power Under Different Models

- Run [`figure_10.R`](code/figure_10.R), [`figure_11.R`](code/figure_11.R), [`figure_12.R`](code/figure_12.R), [`figure_13.R`](code/figure_13.R), [`figure_14.R`](code/figure_14.R) to get [`Figure 10`](figures/figure_power_ratio_delta0.01.pdf), [`Figure 11`](figures/figure_power_ratio_delta0.04.pdf), [`Figure 12`](figures/figure_power_ratio_0.2.pdf), [`Figure 13`](figures/figure_power_ratio_0.3.pdf), [`Figure 14`](figures/figure_power_ratio_0.4.pdf) in the main article, respectively.
  - Contour plots of power ratio between the treatment effect estimator from the HH-ANT and HH models.
 
- Run [`figure_15.R`](code/figure_15.R) and [`figure_16.R`](code/figure_16.R) to get [`Figure 15`](figures/figure_power_ratio_fixed_Delta.pdf) and [`Figure 16`](figures/figure_power_ratio_fixed_ratio.pdf) in the main article, respectively.
  - Contour plots of power ratio between the treatment effect estimator from the ETI-ANT and ETI models.
 
### Reproduce Simulation Studies

- Run [`simulation_I (type I).R`](code/simulation_I%20\(type%20I\).R) and [`simulation_I (power).R`](code/simulation_I%20\(power\).R) to get [`Table 7`](result/simulation_I%20\(type%20I\)/table_7.R) in the main article.

- Run [`simulation_II.R`](code/simulation_II.R) to get [`Table 3`](result/simulation_II/table_3.R) in the main article.

- Run [`simulation_III.R`](code/simulation_III.R) to get [`Table 8`](result/simulation_III/table_8.R) in the main article.

- Run [`simulation_IV.R`](code/simulation_IV.R) to get [`Table 4`](result/simulation_IV/table_4.R) in the main article.

### Trial planning software
We provide an R function [`Trt.Ant.Power`](software/trial_planning.R) to compute statistical power for stepped wedge cluster randomized trials with/without anticipation.

