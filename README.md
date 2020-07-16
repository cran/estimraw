---
output:
  html_document: default
  pdf_document: default
---

# estimraw

<!-- badges: start -->
<!-- badges: end -->

The goal of estimraw is to four-fold table cell frequencies (raw data) from risk ratios (relative risks), risk differences and odds ratios. While raw data can be useful for doing meta-analysis, such data is often not provided by primary studies (with summary statistics being solely presented). Therefore, based on summary statistics (namely, risk ratios, risk differences and odds ratios), this package estimates the value of each cell in a 2x2 table according to the equations described in Di Pietrantonj C (2006) (https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.2287).

The following notation is used (see also table below):

- a: Patients who are exposed/treated and develop the outcome event.
- b: Patients who are exposed/treated and do not develop the outcome event.
- c: Patients who are exposed/treated and do not develop the outcome event.
- d: Patients who are not exposed/treated and do not develop the outcome event.

|             | Event | No event | Total |
|:------------|:-----:|:--------:|:-----:|
| Exposed     |   a   |     b    |  m1   |
| Not exposed |   c   |     d    |  m2   |
| Total       |  e1   |          |       |


Estimating raw data from risk differences and odds ratio involves application of the quadratic formula. Therefore, for such effect sizes, two results sets ("solutions") may be presented, particularly, when there is no information on the total number of participants developing the outcome event (e1). In those cases, the researcher must choose the most adequate solution based on prior clinical knowledge.

The precision with which effect sizes are provided conditions the uncertainty underlying raw data estimation. For example, effect sizes ranging from 1.45 to 1.54 can be reported as '1.5'. Therefore, when e1 is unknown, this function provides not only point estimates based on exact input values, but also a range of compatible values for a, b, c and d (considering possible roundings of effect sizes).  When e1 is known, this function only provides the result(s) in which a+c=e1.


## Installation

You can install the released version of estimraw from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("estimraw")
```

## Getting started

estimraw allows for the use of `estim_raw`function, which has eight arguments:

- `es`: Value of the effect size (summary statistic).
- `lb`: Lower bound of the 95% confidence interval of the effect size.
- `ub`: Upper bound of the 95% confidence interval of the effect size.
- `m1`: Total number of participants in the exposed/treated group.
- `m2`: Total number of participants in the unexposed/control group.
- `e1`: Total number of participants developing the outcome event (optional).
- `dec`: Number of decimal places with which effect sizes are being presented.
- `measure`: Character string indicating the type of effect size measure. Possible options are the risk ratio/relative risk ("rr"), the risk difference ("rd"), and the odds ratio ("or").

The following outputs are obtained:

- Estimates from risk ratios: A dataframe. If there is information on the number of participants developing the outcome event (e1), this dataframe will list all sets of results in which a+c=e1. If no such information is provided, this dataframe will list a point estimate for each cell (calculated based on the exact input values), as well as a minimum and a maximum estimate.
- Estimates from risk differences and odds ratios with known number of participants developing the outcome event (e1): A dataframe containing all sets of results in which a+c=e1.
- Estimates from risk differences and odds ratios with unknown number of participants developing the outcome event (e1): A list consisting of two dataframes - solution1 (presenting the results of the first solution of the quadratic formula) and solution2 (presenting the results of the second solution of the quadratic formula). Each of these dataframes will list a point estimate for each cell (calculated based on the exact input values), as well as a minimum and a maximum estimate.



## Examples

### Example 1 - low effect size precision and no information on the number of participants developing the outcome of event:

Consider a primary study presenting the Risk Ratio as 0.6, with a 95% confidence interval of 0.4-0.9 (worked example 3.1. of Di Petrantonj C (Statistics in Medicine 2006;25:2299-2322)).
Let us assume that (i) We know that there are 352 participants in the treated/exposed group, and 376 participants in the control group; (ii) we do not know the total number of participants developing the outcome of event.

Raw data can be estimated by:

``` r
estim_rr <- estim_raw(es=0.6,lb=0.4,ub=0.9,m1=352,m2=376,dec=1,measure="rr")
estim_rr
```

The results indicate that the value of a (corresponding to the number of treated patients who develop the outcome event) lies between 22 and 47, and that the value of c (corresponding to the number of control patients who develop the outcome event) lies between 38 and 87. This large range reflects the low estimate precision and lack of information on the number of participants developing the outcome event. The real values are actually 37 and 65.


### Example 2 - high effect size precision and no information on the number of participants developing the outcome of event:

Consider a primary study presenting the Odds Ratio as 0.6207, with a 95% confidence interval of 0.3382-1.1391 (worked example 4.2. of Di Pietrantonj C (Statistics in Medicine 2006;25:2299-2322)).
Let us assume that (i) We know that there are 355 participants in the treated/exposed group, and 366 participants in the control group; (ii) we do not know the total number of participants developing the outcome of event.

Raw data can be estimated by:

``` r
estim_or1 <- estim_raw(es=0.6207,lb=0.3382,ub=1.1391,m1=355,m2=366,dec=4,measure="or")
estim_or1
```

The results indicate that the values of a and c are either respectively 18 and 29 (solution 1) or 327 and 348 (solution 2). These two solutions are presented since the quadratic formula had been applied to estimate raw data. The correct solution should be assumed based on prior clinical knowledge (i.e., whether the outcome event is a probable one or not). The real values are actually a=18 and c=29.


### Example 3 - low effect size precision and information on the number of participants developing the outcome of event:

Consider a primary study presenting the Odds Ratio as 0.6, with a 95% confidence interval of 0.3-1.1 (worked example 4.2. of Di Pietrantonj C (Statistics in Medicine 2006;25:2299-2322)). 
Let us assume that (i) We know that there are 355 participants in the treated/exposed group, and 366 participants in the control group; (ii) we know that, overall, 47 participants developed the outcome of event.

Raw data can be estimated by:

``` r
estim_or2 <- estim_raw(es=0.6,lb=0.3,ub=1.1,m1=355,m2=366,e1=47,dec=1,measure="or")
estim_or2
```

The results indicate that the values of a and c are either respectively 17 and 30 or 18 and 29. The real values are actually a=18 and c=29.


### Example 4 - high effect size precision and information on the number of participants developing the outcome of event:

Consider a primary study presenting the risk difference as -7.83%, with a 95% confidence interval of -13.87;-1.8% (worked example 2.2. of Di Pietrantonj C (Statistics in Medicine 2006;25:2299-2322)).
Let us assume that (i) We know that there are 373 participants in the treated/exposed group, and 357 participants in the control group; (ii) we know that, overall, 163 participants developed the outcome of event.

Raw data can be estimated by:

``` r
estim_rd <- estim_raw(es=-0.0783,lb=-0.1387,ub=-0.018,m1=373,m2=357,e1=163,dec=4,measure="rd")
estim_rd
```
The results indicate that the values of a and c are respectively 69 and 94. These actually correspond to the real values of a and c.


## References

- Di Pietrantonj C. Four-fold table cell frequencies imputation in meta analysis. Statistics in Medicine 2006;25:2299-2322.
