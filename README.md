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

## Installation

After the GitHub repository is renamed to `spectropath`, install directly with:

```r
install.packages("remotes")
remotes::install_github("RafaelSdeSouza/spectropath")
```

From a local checkout:

```r
remotes::install_local("/Users/rd23aag/Documents/GitHub/spectropath")
```

During development:

```r
devtools::load_all("/Users/rd23aag/Documents/GitHub/spectropath")
```

## Core path-signature features

```r
library(spectropath)

u <- seq(-5, 5, length.out = 300)
f <- exp(-0.5 * (u / 0.8)^2) + 0.2 * exp(-0.5 * ((u - 2) / 1.2)^2)
path <- cbind(u, f)

path_features(path, depth = 4)
```

`path_features()` returns the compact paper notation by default:

| Column | Scalar helper | Intuition |
|---|---|---|
| `p2` | `levy_area()` | signed area / blue--red handedness |
| `p3u` | `curl_u()` | where asymmetry lies along the ordered coordinate |
| `p3F` | `curl_f()` | whether asymmetry lies in bright or faint parts of the line |
| `p4F` | `jerk_f()` | higher-order flux modulation, shoulders, broad bases |
| `p4T` | `twist_uf()` | twist-like bends, reversals, double peaks, multi-components |
| `p_pm` | `emabs_area()` | emission--absorption ordering |

For readable column names matching the scalar helpers, use:

```r
path_features(path, depth = 4, notation = "descriptive")
```

Additional exploratory contrasts are available with `extended = TRUE`.

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

## Scope

`spectropath` is intentionally narrower than the exploratory `luma` package.
It is meant to make the paper examples reproducible with clean notation and a
stable public API.
