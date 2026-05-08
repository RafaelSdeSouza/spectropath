# spectropath

`spectropath` is a small companion R package for the paper workflow on
spectral lines as ordered paths. It treats a spectral feature as a path

```text
X = (u, F)
```

where `u` is an ordered coordinate, such as wavelength, rest-frame wavelength,
phase, time, or velocity relative to a transition, and `F` is the measured
signal or line flux.

The package exposes only the paper-facing path-signature diagnostics, plus a
small set of classical line-profile summaries for benchmarking and
machine-learning workflows.

## Core path-signature features

```r
library(spectropath)

u <- seq(-5, 5, length.out = 300)
f <- exp(-0.5 * (u / 0.8)^2) + 0.2 * exp(-0.5 * ((u - 2) / 1.2)^2)
path <- cbind(u, f)

path_features(path, depth = 4)
```

The returned columns follow the paper notation:

| Function | Symbol | Intuition |
|---|---:|---|
| `levy_area()` | \(A_{uF}\) | signed area / first order-sensitive asymmetry |
| `curl_u()` | \(C_u\) | where signed area sits along the ordered coordinate |
| `curl_f()` | \(C_F\) | whether signed area sits in core or wings |
| `skew_signature()` | \(S_{\rm skew}\) | order-sensitive blue--red or early--late asymmetry |
| `jerk_u()` | \(J_u\) | higher-order coordinate modulation |
| `twist_uf()` | \(T_{uF}\) | bends, shoulders, reversals, multi-component structure |
| `jerk_f()` | \(J_F\) | higher-order flux modulation |
| `emabs_area()` | \(A_{+-}\) | emission--absorption ordering |

For mixed emission--absorption profiles:

```r
f_mixed <- exp(-0.5 * ((u - 0.7) / 0.5)^2) -
  0.8 * exp(-0.5 * ((u + 0.6) / 0.35)^2)

emabs_area(cbind(u, f_mixed))
```

## Classical convenience features

Classical descriptors are provided for comparison, not as path-signature
quantities:

```r
classical_features(path)
```

This includes integrated line flux, pseudo-equivalent width, peak flux, FWHM,
\(W_{80}\), centroid, line dispersion, skewness, and kurtosis excess.

## Installation from source

From the package directory:

```r
install.packages(".", repos = NULL, type = "source")
```

or during development:

```r
devtools::load_all()
```

## Scope

`spectropath` is intentionally narrower than the exploratory `luma` package.
It is meant to make the paper examples reproducible with clean notation and a
stable public API.
