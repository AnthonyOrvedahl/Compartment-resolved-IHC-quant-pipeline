#!/usr/bin/env python
"""Aggregate run_cohort.py outputs into per-mouse endpoints and statistics.

Deterministic: no randomness, no model calls. Same inputs give the same tables.

Inputs
  --results-dir   directory holding run_cohort.py outputs, one set per slide tag:
                    cells_<tag>.csv  or  cells_<tag>_w<k>.csv   (one per worker)
  --manifest      CSV with columns: slide, sample, group, day
                    slide  must equal the --tag used with run_cohort.py
                    group  e.g. sgEGFP | sgTMPRSS2   (exactly two groups)
                    day    e.g. D3 | D7
  --out-dir       where the tables are written
  --per-mouse     alternative to --results-dir/--manifest: start from an existing
                  per_mouse_all_compartments.csv and recompute only the two
                  statistics tables. Needs no cell-level data.

Rules (each is a constant below, not a judgement made at run time)
  * Compartment of a cell = its QuPath classification, EXCEPT cells classed
    "Epithelial", which are split by distance to the nearest accepted airway:
        d_airway_um <= AIRWAY_BAND_UM  -> "Airway epithelium"
        otherwise                      -> "Alveolar+consolidated epithelium"
  * A cell is DAB-positive when "Cell: DAB OD mean" >= DAB_THRESHOLD. The
    threshold is fixed for the whole cohort, never re-derived per slide.
  * Luminal and indeterminate objects (EXCLUDED) are dropped BEFORE any
    denominator is formed. They are blood, exudate and debris, not cells whose
    infection status means anything.
  * The unit of analysis is the mouse. Cells within a mouse are not independent.
  * No minimum-cell floor is applied. Per-mouse precision is reported as a Wilson
    interval so low-n mice are visible instead of silently dropped.

Outputs (written to --out-dir)
  per_mouse_all_compartments.csv   n, n_pos, pct_pos and Wilson 95% CI per mouse x compartment
  immune_infiltrate_per_sample.csv lymphocytic-infiltrate cells as % of analysed cells, per mouse
  group_contrasts.csv              Mann-Whitney between the two groups, per compartment x stratum
  compartment_contrasts.csv        paired Wilcoxon between compartments, Holm-adjusted
  summary.json                     Friedman test and counts

Statistics notes
  * group_contrasts.csv reports min_possible_p, the smallest two-sided p the test
    can return at those group sizes. When it exceeds 0.05 the test cannot reach
    significance whatever the data, so "n.s." there means no power, not equality.
  * All p-values are unadjusted except p_holm.
"""
import argparse
import glob
import json
import os
from itertools import combinations
from math import comb, sqrt

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, mannwhitneyu, wilcoxon

AIRWAY_BAND_UM = 30.0
DAB_THRESHOLD = 0.19161676646706588
EXCLUDED = ["Airway lumen", "Vessel lumen", "alveolar lumen", "Indeterminate"]
COMPARTMENTS = ["Airway epithelium", "Alveolar+consolidated epithelium",
                "Alveolar_Septa", "Lymphocytic_infiltrate"]
Z95 = 1.959963984540054


def wilson(k, n, z=Z95):
    """Wilson score interval, returned as percentages."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


def min_possible_p(na, nb):
    """Smallest attainable two-sided exact Mann-Whitney p for group sizes na, nb."""
    return 2.0 / comb(na + nb, min(na, nb))


def holm(p):
    """Holm step-down adjusted p-values, in the input order."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, p[i] * (m - k))
        adj[i] = min(run, 1.0)
    return adj


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def load_cells(results_dir, manifest):
    frames = []
    for slide in manifest.slide:
        files = sorted(glob.glob(os.path.join(results_dir, "cells_%s_w*.csv" % slide)))
        single = os.path.join(results_dir, "cells_%s.csv" % slide)
        if not files and os.path.exists(single):
            files = [single]
        if not files:
            raise SystemExit("no cells_%s*.csv in %s (manifest lists slide %s)"
                             % (slide, results_dir, slide))
        d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        d["slide"] = slide
        if d["Object ID"].duplicated().any():
            raise SystemExit("slide %s: %d duplicate Object IDs across worker files - "
                             "a tile was run twice" % (slide, int(d["Object ID"].duplicated().sum())))
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def assign(cells, dab_threshold, band_um):
    c = cells.copy()
    epi = c["Classification"] == "Epithelial"
    c["is_airway_epi"] = epi & (c["d_airway_um"] <= band_um)
    c["compartment"] = np.where(~epi, c["Classification"],
                                np.where(c["is_airway_epi"], "Airway epithelium",
                                         "Alveolar+consolidated epithelium"))
    c["DAB_pos"] = c["Cell: DAB OD mean"] >= dab_threshold
    return c


def per_mouse(analysed, manifest):
    keys = ["slide", "sample", "group", "day"]
    a = analysed.merge(manifest, on="slide", how="left")
    t = (a[a.compartment.isin(COMPARTMENTS)]
         .groupby(keys + ["compartment"])
         .agg(n=("DAB_pos", "size"), n_pos=("DAB_pos", "sum")).reset_index())
    t["pct_pos"] = 100 * t.n_pos / t.n
    t["ci_lo"], t["ci_hi"] = wilson(t.n_pos, t.n)
    t["ci_halfwidth"] = (t.ci_hi - t.ci_lo) / 2
    inf = (a.assign(is_inf=a.Classification == "Lymphocytic_infiltrate")
           .groupby(keys).agg(n_tissue=("is_inf", "size"), n_inf=("is_inf", "sum")).reset_index())
    inf["pct_inf"] = 100 * inf.n_inf / inf.n_tissue
    inf["ci_lo"], inf["ci_hi"] = wilson(inf.n_inf, inf.n_tissue)
    return t, inf


def group_tests(pm, group_a, group_b):
    rows = []
    for comp in COMPARTMENTS:
        s = pm[pm.compartment == comp]
        for label, sub in [(d, s[s.day == d]) for d in sorted(s.day.unique())] + [("pooled", s)]:
            a = sub[sub.group == group_a].pct_pos
            b = sub[sub.group == group_b].pct_pos
            if len(a) < 1 or len(b) < 1:
                continue
            mp = min_possible_p(len(a), len(b))
            rows.append(dict(compartment=comp, stratum=label, n_a=len(a), n_b=len(b),
                             med_a=a.median(), med_b=b.median(),
                             p=mannwhitneyu(a, b).pvalue, min_possible_p=mp,
                             can_reach_005=mp <= 0.05))
    out = pd.DataFrame(rows)
    return out.rename(columns={"n_a": "n_" + group_a, "n_b": "n_" + group_b,
                               "med_a": "med_" + group_a, "med_b": "med_" + group_b})


def compartment_tests(pm):
    w = pm.pivot_table(index="slide", columns="compartment", values="pct_pos")[COMPARTMENTS].dropna()
    rows = []
    for ca, cb in combinations(COMPARTMENTS, 2):
        rows.append(dict(A=ca, B=cb, median_A=w[ca].median(), median_B=w[cb].median(),
                         n_pairs=len(w), p=wilcoxon(w[ca], w[cb]).pvalue))
    t = pd.DataFrame(rows)
    t["p_holm"] = holm(t.p)
    t["stars"] = [stars(p) for p in t.p_holm]
    chi2, pf = friedmanchisquare(*[w[c] for c in COMPARTMENTS])
    return t, dict(n_mice_all_compartments=len(w), friedman_chi2=float(chi2), friedman_p=float(pf))


def round_stats(gt, ct):
    """Human-readable rounding for the two statistics tables (counts stay exact)."""
    if gt is not None:
        gt = gt.copy()
        for col in [c for c in gt.columns if c.startswith("med_")]:
            gt[col] = gt[col].round(2)
        gt["p"] = gt["p"].round(4)
        gt["min_possible_p"] = gt["min_possible_p"].round(4)
    if ct is not None:
        ct = ct.copy()
        ct["median_A"] = ct["median_A"].round(2)
        ct["median_B"] = ct["median_B"].round(2)
        ct["p"] = ct["p"].round(5)
        ct["p_holm"] = ct["p_holm"].round(4)
    return gt, ct


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results-dir")
    ap.add_argument("--manifest")
    ap.add_argument("--per-mouse")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--dab-threshold", type=float, default=DAB_THRESHOLD)
    ap.add_argument("--airway-band-um", type=float, default=AIRWAY_BAND_UM)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    if a.per_mouse:
        pm = pd.read_csv(a.per_mouse, dtype={"slide": str, "sample": str, "group": str, "day": str})
        groups = sorted(pm.group.unique())
        if len(groups) != 2:
            raise SystemExit("per-mouse table must contain exactly two groups, found %s" % groups)
        gt = round_stats(group_tests(pm, groups[0], groups[1]), None)[0]
        ct, summ = compartment_tests(pm)
        ct = round_stats(None, ct)[1]
        gt.to_csv(os.path.join(a.out_dir, "group_contrasts.csv"), index=False)
        ct.to_csv(os.path.join(a.out_dir, "compartment_contrasts.csv"), index=False)
        print("wrote group_contrasts.csv and compartment_contrasts.csv to %s" % a.out_dir)
        return
    if not (a.results_dir and a.manifest):
        raise SystemExit("give --results-dir and --manifest, or --per-mouse")

    man = pd.read_csv(a.manifest, dtype=str)
    need = {"slide", "sample", "group", "day"}
    if not need <= set(man.columns):
        raise SystemExit("manifest needs columns %s, has %s" % (sorted(need), list(man.columns)))
    groups = sorted(man.group.unique())
    if len(groups) != 2:
        raise SystemExit("manifest must contain exactly two groups, found %s" % groups)

    cells = assign(load_cells(a.results_dir, man), a.dab_threshold, a.airway_band_um)
    analysed = cells[~cells.Classification.isin(EXCLUDED)]
    pm, inf = per_mouse(analysed, man)
    gt = group_tests(pm, groups[0], groups[1])
    ct, summ = compartment_tests(pm)

    gt, ct = round_stats(gt, ct)

    pm.to_csv(os.path.join(a.out_dir, "per_mouse_all_compartments.csv"), index=False)
    inf.to_csv(os.path.join(a.out_dir, "immune_infiltrate_per_sample.csv"), index=False)
    gt.to_csv(os.path.join(a.out_dir, "group_contrasts.csv"), index=False)
    ct.to_csv(os.path.join(a.out_dir, "compartment_contrasts.csv"), index=False)
    summ.update(n_cells_detected=int(len(cells)), n_cells_analysed=int(len(analysed)),
                n_excluded_luminal_or_indeterminate=int(len(cells) - len(analysed)),
                groups=groups, dab_threshold=a.dab_threshold, airway_band_um=a.airway_band_um)
    json.dump(summ, open(os.path.join(a.out_dir, "summary.json"), "w"), indent=2)
    print("cells detected %d | analysed %d | excluded %d (%.1f%%)"
          % (len(cells), len(analysed), len(cells) - len(analysed),
             100 * (len(cells) - len(analysed)) / len(cells)))
    print("wrote 4 tables + summary.json to %s" % a.out_dir)


if __name__ == "__main__":
    main()
