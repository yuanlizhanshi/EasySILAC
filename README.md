# EasySILAC: SILAC protein turnover analysis with R and Rust fitting engines

EasySILAC provides the core functions for pulsed-SILAC protein turnover analysis:

- **Quality control** of peptide-level SILAC data ([`quality_control_SILAC`](reference))
- **Steady-state fitting** of heavy-ratio kinetics (`fit_silac`)
- **Non-steady-state fitting** with total-protein dynamics (`fit_silac_nonsteady`)
- **Global labeling-efficiency correction** with a Hill function (`fit_hill_labelling_efficiency`)
- **Visualization** of fitted curves (`plot_protein_silac_fit`, `plot_protein_silac_nonsteady_fit`, `plot_protein_silac_nonsteady_adjust_fit`, ...)

Fitting is available in two numerically identical engines:

- **R** (default): the original `lm`/`nls` implementation
- **Rust** (optional, much faster for large datasets): compiled on first use via
  [rextendr](https://extendr.github.io/rextendr/), requires a
  [Rust toolchain](https://rustup.rs/). Enable globally with
  `silac_use_rust(TRUE)`, or call `fit_silac_rust()` / `fit_silac_batch()` /
  `fit_silac_nonsteady_batch()` directly. The batch variants fit a whole matrix
  multithreaded in a single call.

## Installation

EasySILAC depends on [EasyProtein](https://github.com/yuanlizhanshi/EasyProtein):

``` r
install.packages("devtools")
devtools::install_github("yuanlizhanshi/EasyProtein")
devtools::install_github("yuanlizhanshi/EasySILAC")
```

For the Rust fitting engine additionally install rextendr (and a Rust toolchain):

``` r
install.packages("rextendr")
```

## Quick start

``` r
library(EasySILAC)

# Steady-state fit of one protein's heavy-ratio time course
fit_silac(ratio = c(0.05, 0.12, 0.22, 0.35),
          t = c(1, 2, 4, 6),
          Peptide = "MYGENE")

# Same fit with the Rust engine
fit_silac_rust(ratio = c(0.05, 0.12, 0.22, 0.35),
               t = c(1, 2, 4, 6),
               Peptide = "MYGENE")

# Or switch engines globally for the whole session
silac_use_rust(TRUE)
fit_silac(ratio = c(0.05, 0.12, 0.22, 0.35),
          t = c(1, 2, 4, 6),
          Peptide = "MYGENE")  # now runs in Rust

# Batch fit of a whole ratio matrix (multithreaded Rust)
fit_silac_batch(mtx_ratio, t = t_vec, mtx_total = mtx_total)

# Non-steady-state fit (degradation + synthesis from total intensity)
fit_silac_nonsteady(ratio = ratio_vec, total = total_vec, t = t_vec,
                    Peptide = "MYGENE")

# Global labeling-efficiency (Hill) correction across H100/H75/H50 samples
hill <- fit_hill_labelling_efficiency(global_full_ratio)
```

See `vignette("quickstart", package = "EasySILAC")` for a complete worked example.

## Benchmark

Measured on a 12-core machine with real data (SILAC D0, 160k peptides x 16
samples; 7.7k proteins):

| Task | R engine | Rust engine | Speedup |
|------|----------|-------------|---------|
| Peptide steady-state fit, one concentration (160,625 peptides) | 544 s | 0.26 s | **~2000x** |
| Full peptide fitting step, one day x 3 concentrations | 34.9 min | 76 s | **~28x** |
| Full pipeline (peptide + protein + Hill + adjusted fits), one day | ~45-50 min | 11.4 min | **~4x** |

The Rust engine is numerically identical to the R engine (max difference
~1e-14 on all fit statistics, verified on the full D0 dataset), so results are
interchangeable.

## Note

`geom_mean` and `df2mtx` are re-exported from EasyProtein. EasySILAC keeps its
own `saveplot` (adds `target_dir`) and `plot_pca` (adds `shape_by`), which
extend the EasyProtein versions.
