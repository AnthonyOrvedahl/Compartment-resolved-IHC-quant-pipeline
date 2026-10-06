#!/usr/bin/env python
"""QC: were the exported DAB values computed with the cohort stain vectors?

The positivity threshold is applied to the QuPath measurement "Cell: DAB OD
mean", which depends on the colour-deconvolution vectors set on each image. If one
slide was detected under different vectors (for example QuPath's default H-DAB),
its optical densities are not comparable and the fixed threshold means something
different there. This script recomputes the cell DAB OD from the raw pixels under
each candidate vector set and reports which one reproduces the exported values.

Usage
    python check_stains.py --slide S.ndpi --meas S.Classifications.txt

Output: n cells used, and for each vector set the ratio sum(recomputed)/sum(exported)
and the correlation. The matching set has ratio near 1. Exit status 1 when the
cohort vectors do not match (ratio outside 0.90-1.10).

Method: a fixed random sample (seed 0) of cells with area 30-150 um^2; DAB OD
averaged over a disc of equal area at the exported centroid, level 0; optical
density = -log10(max(value,1)/255); deconvolution with the normalised stain vectors
and their cross-product as the residual. The disc approximates QuPath's cell
outline, so ratios are accurate to a few percent, which is ample to separate
vector sets that differ by 35 percent or more.
"""
import argparse
import sys

import numpy as np
import openslide
import pandas as pd

VECTORS = {
    "cohort": ([0.710096662979583, 0.6422769816906444, 0.2885186441386696],
               [0.436680328623905, 0.5856604407389305, 0.6828705138943986]),
    "qupath_default": ([0.6511078257574492, 0.7011930431234068, 0.29049426072255424],
                       [0.26916687204956063, 0.5682411743268502, 0.777593185920953]),
}
N_CELLS = 400
TOLERANCE = 0.10


def inverse_matrix(h, d):
    h = np.asarray(h, float) / np.linalg.norm(h)
    d = np.asarray(d, float) / np.linalg.norm(d)
    r = np.cross(h, d)
    r /= np.linalg.norm(r)
    return np.linalg.inv(np.vstack([h, d, r]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--slide", required=True)
    ap.add_argument("--meas", required=True)
    a = ap.parse_args()

    sl = openslide.OpenSlide(a.slide)
    mpp = float(sl.properties["openslide.mpp-x"])
    hdr = pd.read_csv(a.meas, sep="\t", nrows=0).columns
    xc = [c for c in hdr if c.startswith("Centroid X")][0]
    yc = [c for c in hdr if c.startswith("Centroid Y")][0]
    ac = [c for c in hdr if c.startswith("Cell: Area")][0]
    q = pd.read_csv(a.meas, sep="\t", usecols=[xc, yc, ac, "Cell: DAB OD mean"])
    q = q[(q[ac] > 30) & (q[ac] < 150)]
    s = q.sample(min(N_CELLS, len(q)), random_state=0)
    inv = {k: inverse_matrix(*v) for k, v in VECTORS.items()}
    obs, pred = [], {k: [] for k in VECTORS}
    for _, r in s.iterrows():
        rad = np.sqrt(r[ac] / np.pi) / mpp
        cx, cy = r[xc] / mpp, r[yc] / mpp
        h = int(np.ceil(rad)) + 1
        x0, y0 = int(cx) - h, int(cy) - h
        if x0 < 0 or y0 < 0:
            continue
        im = np.asarray(sl.read_region((x0, y0), 0, (2 * h + 1, 2 * h + 1)).convert("RGB"),
                        dtype=np.float32)
        yy, xx = np.mgrid[0:2 * h + 1, 0:2 * h + 1]
        disc = (xx - (cx - x0)) ** 2 + (yy - (cy - y0)) ** 2 <= rad ** 2
        if disc.sum() < 10:
            continue
        od = -np.log10(np.maximum(im[disc], 1.0) / 255.0)
        obs.append(r["Cell: DAB OD mean"])
        for k in VECTORS:
            pred[k].append(float((od @ inv[k])[:, 1].mean()))
    obs = np.array(obs)
    print("cells used: %d" % len(obs))
    ratios = {}
    for k in VECTORS:
        p = np.array(pred[k])
        ratios[k] = p.sum() / obs.sum()
        print("  %-15s ratio %.3f   r %.3f" % (k, ratios[k], np.corrcoef(p, obs)[0, 1]))
    ok = abs(ratios["cohort"] - 1.0) <= TOLERANCE
    best = min(ratios, key=lambda k: abs(ratios[k] - 1.0))
    print("RESULT: %s (best match: %s)" % ("cohort vectors match" if ok else "MISMATCH with cohort vectors", best))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
