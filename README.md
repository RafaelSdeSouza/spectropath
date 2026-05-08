# spectropath paper workspace

This repository is organized around the `spectropath` companion package and
the astronomy paper material on spectral lines as paths.

## Active package

- Package: `spectropath`
- Location: `spectropath/`
- Scope: path-signature diagnostics for astronomical spectral-line morphology.

The package exposes the paper-facing quantities:

| Function | Paper symbol |
|---|---:|
| `levy_area()` | \(A_{uF}\) |
| `curl_u()` | \(C_u\) |
| `curl_f()` | \(C_F\) |
| `skew_signature()` | \(S_{\rm skew}\) |
| `jerk_u()` | \(J_u\) |
| `twist_uf()` | \(T_{uF}\) |
| `jerk_f()` | \(J_F\) |
| `emabs_area()` | \(A_{+-}\) |

Classical descriptors such as equivalent width, FWHM, \(W_{80}\), centroid,
line dispersion, skewness, and kurtosis are included as convenience utilities
for comparison and machine-learning workflows. They are not path-signature
quantities.

## Paper material

- Paper workspace: `path_signatures/paper/`
- Manuscript source and assets: `path_signatures/Geometry_astronomy/`
- Generated paper figures: `path_signatures/paper/figures/`

## Legacy package

The exploratory `luma` package has been archived outside this repository at:

```text
/Users/rd23aag/Documents/GitHub/luma_legacy_package_2026-05-08
```

The active package to publish from this repository is `spectropath`.

## Quick check

From the repository root:

```sh
R CMD build spectropath
R CMD check spectropath_0.1.0.tar.gz --no-manual --no-build-vignettes
```

The current package check passes with `Status: OK`.
