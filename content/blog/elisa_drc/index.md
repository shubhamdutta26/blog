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
excerpt: This theme offers built-in Font Awesome icons for organizing your collection
  of social accounts and their links. Use icons to help visitors find you wherever
  you want to be found, and learn how to show or hide them in your site's header,
  footer, homepage, about page, and contact form.
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

Our data is a dose response ELISA experiment for two different antibody reagents. The data can be found here.


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

Wherever you end up wanting to show your social icons, you'll need to start by setting up the links in your site `config.toml` file. Open that up and scroll down to the `[[params.social]]` section. The start of it looks like this:

## The final plot


``` r
theme_set(theme_bw(base_size = 20))
ggplot(summary, aes(x = primary_mab_conc, 
                    y = mean_od, 
                    group = primary_mab_name, 
                    color = primary_mab_name)) +
  geom_point(shape = 21, size = 3, stroke = 1.5) +
  geom_smooth(method = drc::drm, 
              method.args = list(fct = drc::L.4()),
              se = FALSE, linewidth = 1) +
  geom_errorbar(aes(ymin = mean_od - mean_sd,
                    ymax = mean_od + mean_sd),
                width = 0.1) +
  scale_x_log10() +
  scale_color_manual(values = c("#2F9D72", "#2D3047")) +
  labs(x = "Log of antibody concentration (µg/ml)",
       y = expression(OD[450]),
       color = NULL) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.7))
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

<img src="{{< blogdown/postref >}}index_files/figure-html/plot-data-1.png" width="672" />

## Statistical analysis of the fit


``` r
raw_data_no_blanks <- subset(raw_data, primary_mab_name != "blank")
drm_model <- drm(formula = od450~primary_mab_conc, 
                 curveid = primary_mab_name, 
                 data = raw_data_no_blanks, 
                 fct = LL.4(names=c("Slope", "Lower", "Upper", "ED50")))
# tidy(drm_model)
# glance(drm_model)
# augment(drm_model, data = raw_data_no_blanks)
```


``` r
ggplot(augment(drm_model, data = raw_data_no_blanks), 
       aes(.fitted, .resid, color = as.factor(primary_mab_conc))) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 21, size = 3, stroke = 1.5) +
  labs(color = "dilutions") +
  facet_grid(cols = vars(primary_mab_name))
```

```
## Warning in model$der: partial match of 'der' to 'deriv1'
```

<img src="{{< blogdown/postref >}}index_files/figure-html/fit-data-2-1.png" width="672" />
