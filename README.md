A Bayesian Framework for Baseflow Estimation in Rivers
================
Edgar Santos-Fernandez, PhD

This tutorial shows how to estimate baseflow using different methods,
including the Lyne and Hollick (Lyne and Hollick (1979), Ladson et al.
(2013)) and Eckhardt’s filters (Eckhardt (2008)) in R. We use Bayesian
inference in Stan with the `bayesflow` package.

This approach allows:

- incorporating covariates e.g. precipitation

- treat the filter parameters probabilistically and obtain probabilistic
  estimates of the baseflow,

- capture and propagate uncertainty, and better understand the
  underlying processes.

- imputing missing data and fills the gaps on the go.

# Methods of baseflow estimation

Baseflow is the portion of streamflow that is sustained by groundwater
discharge, rather than by surface water runoff or precipitation.

The Lyne and Hollick filter separates the streamflow into two
components: a quick response component that is associated with surface
water runoff and a slow response component that is associated with
groundwater discharge. The method involves applying a low pass filter to
the flow data to separate the two components and then estimating the
baseflow by summing the slow response component over a given period of
time.

The filter is generally defined as follows. Let $y_t$ be the water level
at time $t$. The quick flow response is obtained using:

$$
f_t = \max (0, 
\alpha f_{t-1} + \frac{1+\alpha}{2}(y_{t} - y_{t-1}))   
$$

The $\alpha$ parameter determines the smoothing of the filter. The
baseflow is obtained using:

$$
b_t = y_t - f_t
$$

Eckhardt (2008) suggested an improved method that has been reported to
perform well on different scenarios

$$
b_t = \frac{(1-\omega)\alpha b_{t-1} + (1-\alpha)\omega y_t}{1-\alpha \omega}
$$

where $\omega$ is a baseflow index, usually fixed at 0.8, 0.5 or 0.25
depending on the chracteristics of the stream.

## Installation

And the development version from [GitHub](https://github.com/) with:

``` r
devtools::install_github("EdgarSantos-Fernandez/bayesflow")
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-eckhardt2008comparison" class="csl-entry">

Eckhardt, K. 2008. “A Comparison of Baseflow Indices, Which Were
Calculated with Seven Different Baseflow Separation Methods.” *Journal
of Hydrology* 352 (1-2): 168–73.

</div>

<div id="ref-ladson2013standard" class="csl-entry">

Ladson, Anthony Richard, R Brown, B Neal, and R Nathan. 2013. “A
Standard Approach to Baseflow Separation Using the Lyne and Hollick
Filter.” *Australasian Journal of Water Resources* 17 (1): 25–34.

</div>

<div id="ref-lyne1979stochastic" class="csl-entry">

Lyne, V, and M Hollick. 1979. “Stochastic Time-Variable Rainfall-Runoff
Modelling.” In *Institute of Engineers Australia National Conference*,
79:89–93. 10. Institute of Engineers Australia Barton, Australia.

</div>

</div>
