// Rust port of fit_silac() and fit_silac_nonsteady() from functions.R.
// Replicates the R versions numerically:
//  - lm() is replaced by exact closed-form OLS (identical estimates/SE/r2)
//  - nls(algorithm = "port", lower = 0) for s is replaced by the closed-form
//    bounded least-squares solution (the model is linear in s given k, P0),
//    which equals the converged port estimate to within nls tolerance.
//
// Batch variants fit every row of column-major matrices in parallel with
// scoped threads and drop failed rows, matching purrr::map_dfr over the
// per-item functions.

use extendr_api::prelude::*;

const EPS: f64 = 1e-6;
const MIN_K: f64 = 0.0001;
const MIN_POINTS: usize = 3;

struct Ols {
    intercept: f64,
    slope: f64,
    r2: f64,
    se_intercept: f64,
    se_slope: f64,
}

// Closed-form simple linear regression, numerically equivalent to R's lm(y ~ x).
fn ols(x: &[f64], y: &[f64]) -> Ols {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        sxx += dx * dx;
        sxy += dx * (y[i] - my);
    }
    let slope = sxy / sxx;
    let intercept = my - slope * mx;
    let mut rss = 0.0f64;
    let mut tss = 0.0f64;
    for i in 0..x.len() {
        let r = y[i] - (intercept + slope * x[i]);
        rss += r * r;
        let dy = y[i] - my;
        tss += dy * dy;
    }
    let r2 = 1.0 - rss / tss; // NaN when tss == 0, same as summary.lm
    let sigma2 = rss / (n - 2.0); // n >= 3 guaranteed by caller
    Ols {
        intercept,
        slope,
        r2,
        se_intercept: (sigma2 * (1.0 / n + mx * mx / sxx)).sqrt(),
        se_slope: (sigma2 / sxx).sqrt(),
    }
}

// ---------------------------------------------------------------------------
// Steady-state fit (fit_silac)
// ---------------------------------------------------------------------------

struct SsRow {
    k: f64,
    intercept: f64,
    r2: f64,
    se_k: f64,
    se_intercept: f64,
    p0: f64,
    synthesis_rate: f64,
    h0: f64,
    half_life: f64,
    point_number: i32,
}

fn fit_ss_row(ratio: &[f64], t: &[f64], total: Option<&[f64]>) -> Option<SsRow> {
    let mut rr: Vec<f64> = Vec::new();
    let mut tt: Vec<f64> = Vec::new();
    let mut cc: Vec<f64> = Vec::new();
    for i in 0..ratio.len() {
        let ok = ratio[i].is_finite() && t[i].is_finite() && ratio[i] > 0.0 && ratio[i] < 1.0
            && total.map_or(true, |tv| tv[i].is_finite());
        if ok {
            rr.push(ratio[i]);
            tt.push(t[i]);
            if let Some(tv) = total {
                cc.push(tv[i]);
            }
        }
    }
    if rr.len() < MIN_POINTS {
        return None;
    }
    let mut uniq = tt.clone();
    uniq.sort_by(|a, b| a.partial_cmp(b).unwrap());
    uniq.dedup();
    if uniq.len() < MIN_POINTS {
        return None;
    }
    let y: Vec<f64> = rr
        .iter()
        .map(|r| (1.0 - r.clamp(EPS, 1.0 - EPS)).ln())
        .collect();
    let fit = ols(&tt, &y);
    let k = -fit.slope;
    if !k.is_finite() {
        return None;
    }
    let (p0, synthesis_rate) = if !cc.is_empty() {
        let t0 = tt.iter().cloned().fold(f64::INFINITY, f64::min);
        let vals: Vec<f64> = tt
            .iter()
            .zip(cc.iter())
            .filter(|(ti, _)| **ti == t0)
            .map(|(_, ci)| *ci)
            .collect();
        if vals.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            let p0 = vals.iter().sum::<f64>() / vals.len() as f64;
            (p0, k * p0)
        }
    } else {
        (f64::NAN, f64::NAN)
    };
    Some(SsRow {
        k,
        intercept: fit.intercept,
        r2: fit.r2,
        se_k: fit.se_slope,
        se_intercept: fit.se_intercept,
        p0,
        synthesis_rate,
        h0: 1.0 - fit.intercept.exp(),
        half_life: std::f64::consts::LN_2 / k,
        point_number: rr.len() as i32,
    })
}

/// Per-item drop-in replacement for fit_silac() (returns NULL when the R
/// version returns NULL).
#[extendr]
fn fit_silac_rs(ratio: Vec<f64>, t: Vec<f64>, total: Nullable<Vec<f64>>, peptide: &str) -> Robj {
    let total_vec: Option<Vec<f64>> = match total {
        Nullable::NotNull(v) => Some(v),
        Nullable::Null => None,
    };
    match fit_ss_row(&ratio, &t, total_vec.as_deref()) {
        None => ().into_robj(),
        Some(row) => list!(
            Peptide = peptide,
            point_number = row.point_number,
            intercept = row.intercept,
            k = row.k,
            half_life = row.half_life,
            r2 = row.r2,
            se_intercept = row.se_intercept,
            se_k = row.se_k,
            P0 = row.p0,
            synthesis_rate = row.synthesis_rate,
            H0 = row.h0
        )
        .into_robj(),
    }
}

/// Batch version: fits every row of column-major matrices in parallel.
#[extendr]
fn fit_silac_batch_rs(
    ratio_flat: Vec<f64>,
    n_row: i32,
    n_col: i32,
    t: Vec<f64>,
    total_flat: Vec<f64>,
    peptides: Vec<String>,
) -> Robj {
    let n_row = n_row as usize;
    let n_col = n_col as usize;
    let results = parallel_rows(n_row, n_col, |ratio_idx| {
        let ratio: Vec<f64> = ratio_idx.iter().map(|&i| ratio_flat[i]).collect();
        let total: Vec<f64> = ratio_idx.iter().map(|&i| total_flat[i]).collect();
        fit_ss_row(&ratio, &t, Some(&total)).map(|row| {
            (
                row.k,
                row.intercept,
                row.r2,
                row.se_k,
                row.se_intercept,
                row.p0,
                row.synthesis_rate,
                row.h0,
                row.half_life,
                row.point_number,
            )
        })
    });

    let mut peptide_col: Vec<String> = Vec::new();
    let mut point_number: Vec<i32> = Vec::new();
    let mut intercept: Vec<f64> = Vec::new();
    let mut k: Vec<f64> = Vec::new();
    let mut half_life: Vec<f64> = Vec::new();
    let mut r2: Vec<f64> = Vec::new();
    let mut se_intercept: Vec<f64> = Vec::new();
    let mut se_k: Vec<f64> = Vec::new();
    let mut p0: Vec<f64> = Vec::new();
    let mut synthesis_rate: Vec<f64> = Vec::new();
    let mut h0: Vec<f64> = Vec::new();

    for (r, row) in results.into_iter().enumerate() {
        if let Some((rk, ri, rr2, rsk, rsi, rp0, rsyn, rh0, rhl, rpn)) = row {
            peptide_col.push(peptides[r].clone());
            k.push(rk);
            intercept.push(ri);
            r2.push(rr2);
            se_k.push(rsk);
            se_intercept.push(rsi);
            p0.push(rp0);
            synthesis_rate.push(rsyn);
            h0.push(rh0);
            half_life.push(rhl);
            point_number.push(rpn);
        }
    }

    data_frame!(
        Peptide = peptide_col,
        point_number = point_number,
        intercept = intercept,
        k = k,
        half_life = half_life,
        r2 = r2,
        se_intercept = se_intercept,
        se_k = se_k,
        P0 = p0,
        synthesis_rate = synthesis_rate,
        H0 = h0
    )
    .into_robj()
}

// ---------------------------------------------------------------------------
// Non-steady-state fit (fit_silac_nonsteady)
// ---------------------------------------------------------------------------

struct NsRow {
    point_number_total: i32,
    point_number_old: i32,
    p0: f64,
    l0_fit: f64,
    k_raw: f64,
    k: f64,
    k_is_min_capped: bool,
    half_life: f64,
    synthesis_fit: f64,
    synthesis_norm: f64,
    balance: f64,
    r2_old_decay: f64,
    r2_total_fit: f64,
    rmse_total_fit: f64,
    nrmse_total_fit: f64,
}

fn fit_ns_row(ratio: &[f64], total: &[f64], t: &[f64]) -> Option<NsRow> {
    // filter: finite ratio/total/t, total > 0; old = total * (1 - ratio)
    let mut rr: Vec<f64> = Vec::new();
    let mut tt: Vec<f64> = Vec::new();
    let mut cc: Vec<f64> = Vec::new();
    for i in 0..ratio.len() {
        if ratio[i].is_finite() && total[i].is_finite() && t[i].is_finite() && total[i] > 0.0 {
            rr.push(ratio[i]);
            tt.push(t[i]);
            cc.push(total[i]);
        }
    }
    let point_number_total = rr.len() as i32;

    let mut old_t: Vec<f64> = Vec::new();
    let mut old_log: Vec<f64> = Vec::new();
    for i in 0..rr.len() {
        let old = cc[i] * (1.0 - rr[i]);
        if old.is_finite() && old > 0.0 {
            old_t.push(tt[i]);
            old_log.push(old.ln());
        }
    }
    let point_number_old = old_t.len() as i32;

    let mut uniq = old_t.clone();
    uniq.sort_by(|a, b| a.partial_cmp(b).unwrap());
    uniq.dedup();
    if uniq.len() < MIN_POINTS {
        return None;
    }

    let fit_k = ols(&old_t, &old_log);
    let k_raw = -fit_k.slope;
    let l0_fit = fit_k.intercept.exp();
    if !k_raw.is_finite() {
        return None;
    }
    let k = k_raw.max(MIN_K);

    let t0 = tt.iter().cloned().fold(f64::INFINITY, f64::min);
    let p0_vals: Vec<f64> = tt
        .iter()
        .zip(cc.iter())
        .filter(|(ti, _)| **ti == t0)
        .map(|(_, ci)| *ci)
        .collect();
    if p0_vals.is_empty() {
        return None;
    }
    let p0 = p0_vals.iter().sum::<f64>() / p0_vals.len() as f64;
    if !p0.is_finite() || p0 <= 0.0 {
        return None;
    }

    // df_fit_s: t > t0
    let s_t: Vec<f64> = tt
        .iter()
        .zip(cc.iter())
        .filter(|(ti, _)| **ti > t0)
        .map(|(ti, _)| *ti)
        .collect();
    let s_total: Vec<f64> = tt
        .iter()
        .zip(cc.iter())
        .filter(|(ti, _)| **ti > t0)
        .map(|(_, ci)| *ci)
        .collect();
    if s_t.is_empty() {
        return None;
    }

    // Closed-form bounded LS for s: total = A + s * B
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for i in 0..s_t.len() {
        let tau = s_t[i] - t0;
        let e = (-k * tau).exp();
        let a = p0 * e;
        let b = (1.0 - e) / k;
        num += b * (s_total[i] - a);
        den += b * b;
    }
    let s_fit = if den > 0.0 && num.is_finite() {
        (num / den).max(0.0)
    } else {
        f64::NAN
    };

    let (r2_total_fit, rmse_total_fit, nrmse_total_fit) = if s_fit.is_finite() {
        let n = s_total.len() as f64;
        let mean_obs = s_total.iter().sum::<f64>() / n;
        let mut ss_res = 0.0f64;
        let mut ss_tot = 0.0f64;
        for i in 0..s_t.len() {
            let tau = s_t[i] - t0;
            let e = (-k * tau).exp();
            let pred = p0 * e + (s_fit / k) * (1.0 - e);
            ss_res += (s_total[i] - pred).powi(2);
            ss_tot += (s_total[i] - mean_obs).powi(2);
        }
        let r2 = if ss_tot > 0.0 {
            1.0 - ss_res / ss_tot
        } else {
            f64::NAN
        };
        let rmse = (ss_res / n).sqrt();
        (r2, rmse, rmse / mean_obs)
    } else {
        (f64::NAN, f64::NAN, f64::NAN)
    };

    Some(NsRow {
        point_number_total,
        point_number_old,
        p0,
        l0_fit,
        k_raw,
        k,
        k_is_min_capped: k_raw < MIN_K,
        half_life: std::f64::consts::LN_2 / k,
        synthesis_fit: s_fit,
        synthesis_norm: s_fit / p0,
        balance: s_fit / (k * p0),
        r2_old_decay: fit_k.r2,
        r2_total_fit,
        rmse_total_fit,
        nrmse_total_fit,
    })
}

/// Per-item drop-in replacement for fit_silac_nonsteady().
#[extendr]
fn fit_silac_nonsteady_rs(
    ratio: Vec<f64>,
    total: Vec<f64>,
    t: Vec<f64>,
    peptide: &str,
) -> Robj {
    match fit_ns_row(&ratio, &total, &t) {
        None => ().into_robj(),
        Some(row) => list!(
            Peptide = peptide,
            point_number_total = row.point_number_total,
            point_number_old = row.point_number_old,
            P0 = row.p0,
            L0_fit = row.l0_fit,
            k_raw = row.k_raw,
            k = row.k,
            k_is_min_capped = row.k_is_min_capped,
            half_life = row.half_life,
            synthesis_fit = row.synthesis_fit,
            synthesis_norm = row.synthesis_norm,
            balance = row.balance,
            r2_old_decay = row.r2_old_decay,
            r2_total_fit = row.r2_total_fit,
            rmse_total_fit = row.rmse_total_fit,
            nrmse_total_fit = row.nrmse_total_fit
        )
        .into_robj(),
    }
}

/// Batch version of the non-steady-state fit (rows of column-major matrices).
#[extendr]
fn fit_silac_nonsteady_batch_rs(
    ratio_flat: Vec<f64>,
    n_row: i32,
    n_col: i32,
    total_flat: Vec<f64>,
    t: Vec<f64>,
    peptides: Vec<String>,
) -> Robj {
    let n_row = n_row as usize;
    let n_col = n_col as usize;
    let results = parallel_rows(n_row, n_col, |idx| {
        let ratio: Vec<f64> = idx.iter().map(|&i| ratio_flat[i]).collect();
        let total: Vec<f64> = idx.iter().map(|&i| total_flat[i]).collect();
        fit_ns_row(&ratio, &total, &t)
    });

    let mut peptide_col: Vec<String> = Vec::new();
    let mut point_number_total: Vec<i32> = Vec::new();
    let mut point_number_old: Vec<i32> = Vec::new();
    let mut p0: Vec<f64> = Vec::new();
    let mut l0_fit: Vec<f64> = Vec::new();
    let mut k_raw: Vec<f64> = Vec::new();
    let mut k: Vec<f64> = Vec::new();
    let mut k_is_min_capped: Vec<bool> = Vec::new();
    let mut half_life: Vec<f64> = Vec::new();
    let mut synthesis_fit: Vec<f64> = Vec::new();
    let mut synthesis_norm: Vec<f64> = Vec::new();
    let mut balance: Vec<f64> = Vec::new();
    let mut r2_old_decay: Vec<f64> = Vec::new();
    let mut r2_total_fit: Vec<f64> = Vec::new();
    let mut rmse_total_fit: Vec<f64> = Vec::new();
    let mut nrmse_total_fit: Vec<f64> = Vec::new();

    for (r, row) in results.into_iter().enumerate() {
        if let Some(row) = row {
            peptide_col.push(peptides[r].clone());
            point_number_total.push(row.point_number_total);
            point_number_old.push(row.point_number_old);
            p0.push(row.p0);
            l0_fit.push(row.l0_fit);
            k_raw.push(row.k_raw);
            k.push(row.k);
            k_is_min_capped.push(row.k_is_min_capped);
            half_life.push(row.half_life);
            synthesis_fit.push(row.synthesis_fit);
            synthesis_norm.push(row.synthesis_norm);
            balance.push(row.balance);
            r2_old_decay.push(row.r2_old_decay);
            r2_total_fit.push(row.r2_total_fit);
            rmse_total_fit.push(row.rmse_total_fit);
            nrmse_total_fit.push(row.nrmse_total_fit);
        }
    }

    data_frame!(
        Peptide = peptide_col,
        point_number_total = point_number_total,
        point_number_old = point_number_old,
        P0 = p0,
        L0_fit = l0_fit,
        k_raw = k_raw,
        k = k,
        k_is_min_capped = k_is_min_capped,
        half_life = half_life,
        synthesis_fit = synthesis_fit,
        synthesis_norm = synthesis_norm,
        balance = balance,
        r2_old_decay = r2_old_decay,
        r2_total_fit = r2_total_fit,
        rmse_total_fit = rmse_total_fit,
        nrmse_total_fit = nrmse_total_fit
    )
    .into_robj()
}

// ---------------------------------------------------------------------------
// Shared parallel row-map over column-major matrices (extendr types are not
// Send, so worker threads return plain Rust structs and the R object is built
// on the main thread).
// ---------------------------------------------------------------------------
fn parallel_rows<T, F>(n_row: usize, n_col: usize, f: F) -> Vec<Option<T>>
where
    T: Send,
    F: Fn(&[usize]) -> Option<T> + Sync,
{
    let n_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(12);
    let chunk = n_row.div_ceil(n_threads);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..n_row)
            .step_by(chunk)
            .map(|start| {
                let end = (start + chunk).min(n_row);
                let f = &f;
                scope.spawn(move || {
                    (start..end)
                        .map(|r| {
                            let idx: Vec<usize> = (0..n_col).map(|c| r + c * n_row).collect();
                            f(&idx)
                        })
                        .collect::<Vec<Option<T>>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap())
            .collect::<Vec<Vec<Option<T>>>>()
            .into_iter()
            .flatten()
            .collect()
    })
}

extendr_module! {
    mod rextendr;
    fn fit_silac_rs;
    fn fit_silac_nonsteady_rs;
    fn fit_silac_batch_rs;
    fn fit_silac_nonsteady_batch_rs;
}
