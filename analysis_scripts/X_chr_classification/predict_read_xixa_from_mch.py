#!/usr/bin/env python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln, logsumexp


def logit_fraction(values, n_sites):
    adjusted = (values * n_sites + 0.5) / (n_sites + 1.0)
    return np.log(adjusted / (1.0 - adjusted))


def beta_binomial_logpmf(k, n, alpha, beta):
    return (
        gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
        + betaln(k + alpha, n - k + beta) - betaln(alpha, beta)
    )


def initial_beta_parameters(fractions):
    fractions = np.asarray(fractions, dtype=float)
    mean = np.clip(np.mean(fractions), 1e-4, 1 - 1e-4)
    variance = max(np.var(fractions), 1e-5)
    concentration = np.clip(mean * (1 - mean) / variance - 1, 2.0, 500.0)
    return mean * concentration, (1 - mean) * concentration


def fit_beta_binomial_mixture(total_mch, n_ch, n_components):
    total_mch = np.asarray(total_mch, dtype=float)
    n_ch = np.asarray(n_ch, dtype=float)
    fractions = total_mch / n_ch

    def unpack(theta):
        alpha = np.exp(theta[:n_components])
        beta = np.exp(theta[n_components:2 * n_components])
        if n_components == 1:
            weights = np.ones(1)
        else:
            first_weight = expit(theta[-1])
            weights = np.array([first_weight, 1 - first_weight])
        return alpha, beta, weights

    def objective(theta):
        alpha, beta, weights = unpack(theta)
        component_log_prob = np.column_stack([
            np.log(weights[j])
            + beta_binomial_logpmf(total_mch, n_ch, alpha[j], beta[j])
            for j in range(n_components)
        ])
        return -logsumexp(component_log_prob, axis=1).sum()

    starts = []
    if n_components == 1:
        alpha, beta = initial_beta_parameters(fractions)
        starts.append(np.log([alpha, beta]))
    else:
        for split_quantile in (0.4, 0.5, 0.6):
            split = np.quantile(fractions, split_quantile)
            low = fractions <= split
            if low.all() or (~low).all():
                continue
            alpha_low, beta_low = initial_beta_parameters(fractions[low])
            alpha_high, beta_high = initial_beta_parameters(fractions[~low])
            weight_logit = np.log(low.mean() / (1 - low.mean()))
            starts.append(np.array([
                np.log(alpha_low), np.log(alpha_high),
                np.log(beta_low), np.log(beta_high),
                weight_logit,
            ]))
        if not starts:
            raise ValueError("Cannot initialize two beta-binomial components")

    bounds = [(-8, 12)] * (2 * n_components)
    if n_components == 2:
        bounds.append((-8, 8))
    fits = [
        minimize(objective, start, method="L-BFGS-B", bounds=bounds)
        for start in starts
    ]
    fit = min(fits, key=lambda result: result.fun)
    if not np.isfinite(fit.fun):
        raise RuntimeError("Beta-binomial mixture optimization failed")

    alpha, beta, weights = unpack(fit.x)
    component_log_prob = np.column_stack([
        np.log(weights[j])
        + beta_binomial_logpmf(total_mch, n_ch, alpha[j], beta[j])
        for j in range(n_components)
    ])
    log_normalizer = logsumexp(component_log_prob, axis=1)
    posterior = np.exp(component_log_prob - log_normalizer[:, None])
    n_parameters = 2 * n_components + (n_components - 1)
    bic = n_parameters * np.log(total_mch.size) + 2 * fit.fun
    return {
        "alpha": alpha,
        "beta": beta,
        "weights": weights,
        "means": alpha / (alpha + beta),
        "posterior": posterior,
        "bic": bic,
        "converged": fit.success,
    }


def add_figure_columns(df):
    df = df.copy()
    df["span_bp"] = df["max_ref_position"] - df["min_ref_position"] + 1
    df["midpoint"] = ((df["min_ref_position"] + df["max_ref_position"]) / 2).astype(int)
    df["total_mCH_level"] = df["mCH_fraction"]
    df["hmCG_fraction"] = df["n_hmCG"] / df["n_CG_pass"].replace(0, np.nan)
    df["mCG_only_fraction"] = df["n_mCG"] / df["n_CG_pass"].replace(0, np.nan)
    df["hmCG_share_of_mCG"] = df["n_hmCG"] / df["total_mCG"].replace(0, np.nan)
    df["hmCH_fraction"] = df["n_hmCH"] / df["n_CH_pass"].replace(0, np.nan)
    df["mCH_only_fraction"] = df["n_mCH"] / df["n_CH_pass"].replace(0, np.nan)
    return df


def classify_reads(df, args):
    out = add_figure_columns(df)
    out["sample"] = args.sample
    out["sex"] = args.sex
    out["celltype"] = args.celltype
    out["passes_ch_filter"] = out["n_CH_pass"] >= args.min_ch_sites
    out["passes_span_filter"] = out["span_bp"] >= args.min_span_bp
    out["passes_classifier_filter"] = out["passes_ch_filter"] & out["passes_span_filter"]
    out["mCH_class"] = np.where(
        out["passes_classifier_filter"],
        "unclassified_no_bin_fit",
        np.where(out["passes_ch_filter"], "unclassified_short_read", "unclassified_low_CH_sites"),
    )
    out["x_activity_prediction"] = out["mCH_class"]
    out["posterior_probability"] = np.nan
    out["prediction_confidence"] = "unclassified"
    out["classifier_bin_start"] = np.nan
    out["classifier_bin_end"] = np.nan
    out["classifier_status"] = np.where(
        out["passes_classifier_filter"],
        "unclassified_no_bin_fit",
        np.where(out["passes_ch_filter"], "unclassified_short_read", "unclassified_low_CH_sites"),
    )
    out["bb_component"] = np.nan
    out["bb_label"] = pd.NA
    out["bb_posterior"] = np.nan
    out["logit_mCH_fraction"] = logit_fraction(
        out["mCH_fraction"].to_numpy(), out["n_CH_pass"].to_numpy()
    )

    used = out.loc[out["passes_classifier_filter"]].copy()
    if used.empty:
        return out, pd.DataFrame()

    used["classifier_bin_start"] = (used["midpoint"] // args.bin_size) * args.bin_size
    used["classifier_bin_end"] = used["classifier_bin_start"] + args.bin_size
    fit_rows = []

    for bin_start, part in used.groupby("classifier_bin_start", sort=True):
        idx = part.index
        base_row = {
            "sample": args.sample,
            "sex": args.sex,
            "celltype": args.celltype,
            "bin_start": int(bin_start),
            "bin_end": int(bin_start + args.bin_size),
            "n_reads": part.shape[0],
            "min_ch_sites": args.min_ch_sites,
            "min_span_bp": args.min_span_bp,
        }
        if part.shape[0] < args.min_reads:
            used.loc[idx, "classifier_status"] = "unclassified_not_enough_bin_reads"
            fit_rows.append({**base_row, "status": "not_enough_bin_reads"})
            continue

        total_mch = (part["n_mCH"] + part["n_hmCH"]).to_numpy()
        n_ch = part["n_CH_pass"].to_numpy()
        try:
            one = fit_beta_binomial_mixture(total_mch, n_ch, n_components=1)
            two = fit_beta_binomial_mixture(total_mch, n_ch, n_components=2)
        except (RuntimeError, ValueError) as error:
            used.loc[idx, "classifier_status"] = "unclassified_fit_failed"
            fit_rows.append({**base_row, "status": "fit_failed", "fit_error": str(error)})
            continue

        posterior = two["posterior"]
        component = posterior.argmax(axis=1)
        means = two["means"]
        low_component = int(np.argmin(means))
        high_component = int(np.argmax(means))
        low_count = int((component == low_component).sum())
        high_count = int((component == high_component).sum())
        low_fraction = low_count / part.shape[0]
        high_fraction = high_count / part.shape[0]
        minor_fraction = min(low_fraction, high_fraction)
        delta_bic = one["bic"] - two["bic"]
        has_two_components = (
            delta_bic > args.bic_delta
            and minor_fraction >= args.min_component_fraction
        )

        component_label = np.where(component == high_component, "Xa_like", "Xi_like")
        max_post = posterior.max(axis=1)
        used.loc[idx, "bb_component"] = component
        used.loc[idx, "bb_label"] = component_label
        used.loc[idx, "bb_posterior"] = max_post

        fit_rows.append({
            **base_row,
            "status": "fit" if has_two_components else "not_two_component",
            "bic_one_component": one["bic"],
            "bic_two_component": two["bic"],
            "delta_bic_one_minus_two": delta_bic,
            "low_component_mean": means[low_component],
            "high_component_mean": means[high_component],
            "low_component_alpha": two["alpha"][low_component],
            "low_component_beta": two["beta"][low_component],
            "high_component_alpha": two["alpha"][high_component],
            "high_component_beta": two["beta"][high_component],
            "low_component_reads": low_count,
            "high_component_reads": high_count,
            "xa_like_component_fraction": high_fraction,
            "xi_like_component_fraction": low_fraction,
            "minor_component_fraction": minor_fraction,
        })

        if not has_two_components:
            used.loc[idx, "classifier_status"] = "unclassified_not_two_component"
            continue

        confident = max_post >= args.confidence
        used.loc[idx, "mCH_class"] = np.where(
            confident,
            np.where(component == high_component, "high_mCH", "low_mCH"),
            "ambiguous",
        )
        used.loc[idx, "x_activity_prediction"] = np.where(
            confident, component_label, "ambiguous"
        )
        used.loc[idx, "posterior_probability"] = max_post
        used.loc[idx, "prediction_confidence"] = np.where(
            confident, "high", "ambiguous"
        )
        used.loc[idx, "classifier_status"] = "fit"

    assign_cols = [
        "mCH_class", "x_activity_prediction", "posterior_probability",
        "prediction_confidence", "classifier_bin_start", "classifier_bin_end",
        "classifier_status", "bb_component", "bb_label", "bb_posterior",
    ]
    out.loc[used.index, assign_cols] = used[assign_cols]
    return out, pd.DataFrame(fit_rows)


def make_summary(predictions, fit_df, args):
    status = fit_df["status"] if "status" in fit_df else pd.Series(dtype=object)
    rows = []
    for label in ["Xa_like", "Xi_like", "ambiguous"]:
        part = predictions.loc[predictions["x_activity_prediction"].eq(label)]
        rows.append({
            "sample": args.sample,
            "sex": args.sex,
            "celltype": args.celltype,
            "status": "complete",
            "class": label,
            "n_reads": part.shape[0],
            "median_total_mCH_level": part["total_mCH_level"].median(),
            "mean_total_mCH_level": part["total_mCH_level"].mean(),
            "median_CH_sites": part["n_CH_pass"].median(),
            "min_ch_sites": args.min_ch_sites,
            "min_span_bp": args.min_span_bp,
            "bin_size": args.bin_size,
            "min_reads": args.min_reads,
            "min_component_fraction": args.min_component_fraction,
            "confidence_threshold": args.confidence,
            "bic_delta": args.bic_delta,
            "n_bins_fit": int(status.eq("fit").sum()),
            "n_bins_not_two_component": int(status.eq("not_two_component").sum()),
            "n_bins_not_enough_reads": int(status.eq("not_enough_bin_reads").sum()),
            "n_bins_fit_failed": int(status.eq("fit_failed").sum()),
            "assumption": "within_bin_low_total_mCH_component_is_Xi_like",
        })
    return pd.DataFrame(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run per-bin beta-binomial Xi-like/Xa-like read classification."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bin-fits", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sex", required=True)
    parser.add_argument("--celltype", required=True)
    parser.add_argument("--min-ch-sites", type=int, default=100)
    parser.add_argument("--min-span-bp", type=int, default=6000)
    parser.add_argument("--min-reads", type=int, default=100)
    parser.add_argument("--confidence", type=float, default=0.8)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-component-fraction", type=float, default=0.30)
    parser.add_argument("--bic-delta", type=float, default=0.0)
    return parser.parse_args()


def main():
    args = parse_args()
    for path in [args.output, args.bin_fits, args.summary]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    predictions, fit_df = classify_reads(df, args)
    summary = make_summary(predictions, fit_df, args)

    predictions.to_csv(args.output, sep="\t", index=False, compression="infer")
    fit_df.to_csv(args.bin_fits, sep="\t", index=False)
    summary.to_csv(args.summary, sep="\t", index=False)

    print(f"Wrote {args.output}", flush=True)
    print(f"Wrote {args.bin_fits}", flush=True)
    print(f"Wrote {args.summary}", flush=True)


if __name__ == "__main__":
    main()
