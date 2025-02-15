---
author: Shubham Dutta
categories:
- ELISA
- dose response
- ggplot2
- drc
- broom
date: "2023-06-08"
draft: false
excerpt: A step-by-step guide to fitting dose-response curves on ELISA data in R using `ggplot2`, `drc`, and `broom.` This tutorial walks through data preparation, visualization, curve fitting, and statistical analysis of the model.
layout: single
subtitle: A working example of how to fit dose response curves on ELISA data in R.
title: Dose response curve fitting in R
---

Here is an example of how to fit and analyse dose response data using `ggplot2`, `drc`, and `broom`. We will start by loading some packages needed for the analysis.


``` r
library(readr)
library(dplyr)
library(ggplot2)
library(drc)
library(broom)
```

## The data

Our data is a dose response ELISA experiment for two different antibody reagents. The data can be found [here](elisa_drc.csv).


``` r
raw_data <- read_csv("elisa_drc.csv", show_col_types = FALSE)
glimpse(raw_data)
```

```
## Rows: 50
## Columns: 10
## $ well                <chr> "A01", "A02", "A03", "A04", "A05", "A06", "B01", "…
## $ coat_protein_name   <chr> "sBACE", "sBACE", "sBACE", "sBACE", "sBACE", "sBAC…
## $ coat_protein_ug     <dbl> 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, …
## $ coat_protein_source <chr> "NS0", "NS0", "NS0", "NS0", "NS0", "NS0", "NS0", "…
## $ primary_mab_name    <chr> "aBACE", "aBACE", "aBACE", "aSARS", "aSARS", "aSAR…
## $ primary_mab_clone   <chr> "6626.1", "6626.1", "6626.1", "CC6.29", "CC6.29", …
## $ primary_mab_conc    <dbl> 10.0000000, 10.0000000, 10.0000000, 10.0000000, 10…
## $ secondary_mab_name  <chr> "goat-aHuman", "goat-aHuman", "goat-aHuman", "goat…
## $ secondary_mab_dil   <chr> "1to5000", "1to5000", "1to5000", "1to5000", "1to50…
## $ od450               <dbl> 1.396, 1.170, 1.299, 1.324, 1.170, 1.299, 1.374, 1…
```

Let's focus on three variables.
-   `primary_mab_name`: Two antibodies used in the ELISA experiment.
-   `primary_mab_conc`: Antibody concentration in µg/ml (the dose).
-   `od450`: Absorbance at 450 nm (the response).

## Prepare data before plotting

Before plotting, we need to preprocess the data by subtracting the mean blank OD value from each measurement. This ensures that our response values are baseline-corrected, allowing for a more accurate dose-response analysis. The data is then grouped by antibody and concentration to compute the mean and standard deviation for each condition.


``` r
blank_data <- raw_data |>
    dplyr::filter(primary_mab_name == "blank")
mean_blank <- mean(blank_data[["od450"]], na.rm = TRUE)

summary <- raw_data |>
    dplyr::filter(primary_mab_name != "blank") |>
    dplyr::mutate(blanked_od = od450 - mean_blank) |>
    dplyr::group_by(primary_mab_name, primary_mab_conc) |>
    dplyr::summarise(
      mean_od = mean(blanked_od, na.rm = TRUE),
      mean_sd = sd(blanked_od, na.rm = TRUE),
      .groups = 'drop'
    )
head(summary)
```

```
## # A tibble: 6 × 4
##   primary_mab_name primary_mab_conc mean_od mean_sd
##   <chr>                       <dbl>   <dbl>   <dbl>
## 1 aBACE                     0.00244   0.138  0.0797
## 2 aBACE                     0.00977   0.214  0.0626
## 3 aBACE                     0.0391    0.395  0.0644
## 4 aBACE                     0.156     0.740  0.0531
## 5 aBACE                     0.625     1.11   0.0358
## 6 aBACE                     2.5       1.28   0.0486
```

## The final plot

The plot visualizes the dose-response relationship of two antibody reagents using four-parameter logistic (4PL) in the `drc` package.


``` r
theme_set(theme_bw(base_size = 20))
ggplot(summary, aes(x = primary_mab_conc, 
                    y = mean_od, 
                    group = primary_mab_name, 
                    color = primary_mab_name,
                    shape = primary_mab_name)) +
  geom_point(size = 3, stroke = 1.5) +
  geom_smooth(method = drc::drm, 
              method.args = list(fct = drc::L.4()),
              se = FALSE, linewidth = 1) +
  geom_errorbar(aes(ymin = mean_od - mean_sd,
                    ymax = mean_od + mean_sd),
                width = 0.1) +
  scale_x_log10() +
  scale_color_manual(name = NULL, values = c("#2F9D72", "#2D3047")) +
  scale_shape_manual(name = NULL, values = c(1, 2)) +
  labs(x = "Log concentration (µg/ml)",
       y = expression(OD[450]),
       color = NULL) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.7))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/plot-data-1.png" width="672" />

## Statistical analysis of the fit

After fitting the logistic regression model, we analyze the goodness of fit and residuals. The `tidy()`, `glance()`, and `augment()` functions from the `broom` package allow us to:

-   Extract parameter estimates, including the slope and effective dose (ED50).
-   Assess model fit statistics.
-   Examine residuals to check for systematic deviations or patterns.


``` r
raw_data_no_blanks <- subset(raw_data, primary_mab_name != "blank")
drm_model <- drm(formula = od450~primary_mab_conc, 
                 curveid = primary_mab_name, 
                 data = raw_data_no_blanks, 
                 fct = LL.4(names=c("Slope", "Lower", "Upper", "ED50")))
tidy(drm_model)
```

```
## # A tibble: 8 × 6
##   term  curve estimate std.error statistic  p.value
##   <chr> <chr>    <dbl>     <dbl>     <dbl>    <dbl>
## 1 Slope aBACE   -1.06     0.140      -7.56 8.83e- 9
## 2 Slope aSARS   -1.12     0.147      -7.61 7.66e- 9
## 3 Lower aBACE    0.179    0.0376      4.77 3.44e- 5
## 4 Lower aSARS    0.168    0.0355      4.74 3.70e- 5
## 5 Upper aBACE    1.33     0.0327     40.7  9.60e-31
## 6 Upper aSARS    1.31     0.0314     41.7  4.42e-31
## 7 ED50  aBACE    0.134    0.0180      7.46 1.19e- 8
## 8 ED50  aSARS    0.137    0.0175      7.83 4.05e- 9
```

``` r
glance(drm_model)
```

```
## # A tibble: 1 × 4
##     AIC   BIC logLik   df.residual
##   <dbl> <dbl> <logLik>       <int>
## 1 -106. -90.0 61.82541          34
```
