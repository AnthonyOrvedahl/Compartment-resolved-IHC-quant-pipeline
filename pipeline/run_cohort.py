"""Per-slide driver: classify pixels, find airways, give every cell its distance
to the nearest accepted airway. One slide per invocation; shardable by worker.

Usage
    python run_cohort.py --slide S.ndpi --meas S.Classifications.txt \
        --classifier "Lung feature prediction.json" --outdir out --tag 001 \
        [--worker 0 --nworkers 8]

Inputs
    --slide        OpenSlide-readable whole-slide image
    --meas         QuPath detection export, tab-separated. Columns used:
                   Object ID, Classification, Centroid X µm, Centroid Y µm,
                   Cell: DAB OD mean
    --classifier   QuPath OpenCV pixel-classifier JSON (output type PROBABILITY)

Outputs (in --outdir; suffix _w<k> when --nworkers > 1)
    cells_<tag>.csv      Object ID, Classification, Cell: DAB OD mean,
                         d_airway_um, slide
    cavities_<tag>.csv   one row per candidate cavity: area, margin, lymphoid and
                         adipose prominence, call

Method, per tile (TILE_UM square plus a MARGIN_UM border read on every side)
    1. Gaussian (sigma 1 px) on R, G, B at level 0; classifier applied every
       STRIDE-th pixel; class map then registered to the level-1 grid (every
       second level-0 pixel) so it overlays the tissue mask.
    2. Tissue mask at level 1: hematoxylin OD + DAB OD > 0.15 (cohort stain
       vectors below).
    3. airway_acceptance: class prominences, cavity enumeration, acceptance.
    4. Each exported cell in the tile: distance from its centroid to the nearest
       accepted-airway cavity.

Tiles are laid over the bounding box of the exported cells, not the slide, so a
second tissue section on the same slide that the annotation excluded is never
processed. Tiles are independent: a structure crossing a tile edge is scored once
per side. Each cell falls in exactly one tile, so cell assignment is unaffected.

Deterministic: no random number generation and no model calls; the classifier is
a fixed trained forest applied as-is.
"""
import os, sys, json, time, argparse, tempfile
import numpy as np, pandas as pd, cv2, openslide
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import airway_acceptance as aa

TILE_UM, MARGIN_UM, STRIDE = 1400.0, 400.0, 3
AIRWAY_BAND_UM = 30.0

HEMA = np.array([0.71010, 0.64228, 0.28852]); DAB = np.array([0.43668, 0.58566, 0.68287])
_r = np.cross(HEMA, DAB); _r /= np.linalg.norm(_r)
MINV = np.linalg.inv(np.vstack([HEMA, DAB, _r]).astype(np.float32))


def deconv(rgb):
    od = -np.log10((rgb.astype(np.float32) + 1.0) / 256.0)
    return (od.reshape(-1, 3) @ MINV).reshape(rgb.shape)


def load_forest(path):
    j = json.load(open(path))
    classes = [c["name"] for c in j["metadata"]["outputChannels"]]

    def find(n):
        if isinstance(n, dict):
            if "opencv_ml_rtrees" in n:
                return n
            if n.get("type") == "op.ml.opencv-statmodel" and "model" in n:
                return n["model"]["statmodel"]
            for v in n.values():
                r = find(v)
                if r is not None:
                    return r
        elif isinstance(n, list):
            for v in n:
                r = find(v)
                if r is not None:
                    return r
    tmp = tempfile.mktemp(suffix=".json")
    json.dump(find(j["op"]), open(tmp, "w"))
    rt = cv2.ml.RTrees_load(tmp)
    os.unlink(tmp)
    return rt, classes, j["metadata"]["inputResolution"]["pixelWidth"]["value"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slide", required=True)
    ap.add_argument("--classifier", required=True)
    ap.add_argument("--allow-resolution-mismatch", action="store_true")
    ap.add_argument("--meas", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--nworkers", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    t0 = time.time()
    sl = openslide.OpenSlide(a.slide)
    mpp = float(sl.properties["openslide.mpp-x"]); mpp1 = mpp * 2
    rt, classes, clf_mpp = load_forest(a.classifier)
    if abs(mpp / clf_mpp - 1.0) > 0.05 and not a.allow_resolution_mismatch:
        raise SystemExit("slide pixel size %.4f um differs from the classifier's training "
                         "resolution %.4f um by more than 5%%; the classifier would be applied at the "
                         "wrong scale. Re-train at this resolution, or pass "
                         "--allow-resolution-mismatch to override." % (mpp, clf_mpp))

    q = pd.read_csv(a.meas, sep="\t")
    q["_x"] = (q["Centroid X µm"] / mpp).astype(int)
    q["_y"] = (q["Centroid Y µm"] / mpp).astype(int)
    step0 = int(round(TILE_UM / mpp)); mg0 = int(round(MARGIN_UM / mpp))
    bx0, bx1 = int(q._x.min()), int(q._x.max())
    by0, by1 = int(q._y.min()), int(q._y.max())
    tiles = [(tx, ty)
             for ty in range(by0, by1 + 1, step0)
             for tx in range(bx0, bx1 + 1, step0)]
    mine = tiles[a.worker::a.nworkers]
    print("[%s w%d] %d of %d tiles | cells %d | bbox %d-%d x %d-%d px"
          % (a.tag, a.worker, len(mine), len(tiles), len(q), bx0, bx1, by0, by1), flush=True)

    cav_rows, cell_frames, done = [], [], 0
    for tx, ty in mine:
        side = step0 + 2 * mg0
        im = np.asarray(sl.read_region((tx - mg0, ty - mg0), 0, (side, side)).convert("RGB"))
        if im[..., 0].std() < 1.0:
            done += 1
            continue
        feat = np.stack([gaussian_filter(im[..., c].astype(np.float32), 1.0) for c in range(3)], -1)
        sub = feat[::STRIDE, ::STRIDE].reshape(-1, 3)
        votes = np.zeros((len(sub), len(classes)), np.float32)
        for i in range(0, len(sub), 800_000):
            blk = sub[i:i + 800_000]
            cl = rt.predict(blk)[1].ravel().astype(np.int32)
            votes[np.arange(i, i + len(blk)), cl] = 1.0
        hs = feat[::STRIDE, ::STRIDE].shape[:2]
        pred = votes.argmax(1).reshape(hs).astype(np.uint8)
        im1 = im[::2, ::2]
        # pred is on a STRIDE-spaced level-0 grid. Repeat it back to level-0
        # resolution, THEN take every second pixel to land on the level-1 grid of
        # im1/tis1. (Cropping the level-0-resolution map to im1's shape instead
        # covers only the top-left quarter of the tile at twice the scale.)
        pred0 = np.repeat(np.repeat(pred, STRIDE, 0), STRIDE, 1)
        pred1 = pred0[::2, ::2][:im1.shape[0], :im1.shape[1]]
        C = deconv(im1)
        tis1 = (C[..., 0] + C[..., 1]) > 0.15

        proms = aa.class_prominences(pred1, classes, tis1, mpp1)
        lab, ids, sizes = aa.enumerate_cavities(tis1, mpp1)
        rows = aa.score_cavities(lab, ids, sizes, proms, tis1, mpp1)
        core = slice(int(mg0 / 2), int(mg0 / 2) + int(step0 / 2))
        for r in rows:
            cy, cx = ndi.center_of_mass(lab == r["cavity"])
            if core.start <= cy < core.stop and core.start <= cx < core.stop:
                r.update(slide=a.tag, tile_x=tx, tile_y=ty,
                         level0_x=tx - mg0 + int(cx) * 2, level0_y=ty - mg0 + int(cy) * 2)
                cav_rows.append(r)

        m = (q._x >= tx) & (q._x < tx + step0) & (q._y >= ty) & (q._y < ty + step0)
        s = q[m]
        if len(s):
            lx = ((s._x - (tx - mg0)) / 2).clip(0, tis1.shape[1] - 1).astype(int)
            ly = ((s._y - (ty - mg0)) / 2).clip(0, tis1.shape[0] - 1).astype(int)
            aw = [r["cavity"] for r in rows if r["call"] == "airway"]
            awm = np.isin(lab, aw) if aw else np.zeros(lab.shape, bool)
            d_aw = ndi.distance_transform_edt(~awm) * mpp1
            out = s[["Object ID", "Classification", "Cell: DAB OD mean"]].copy()
            out["d_airway_um"] = d_aw[ly.to_numpy(), lx.to_numpy()]
            out["slide"] = a.tag
            cell_frames.append(out)
        done += 1
        if done % 5 == 0:
            print("  [%s w%d] %d/%d tiles, %.0f s, %d cavities"
                  % (a.tag, a.worker, done, len(mine), time.time() - t0, len(cav_rows)), flush=True)

    sfx = "_w%d" % a.worker if a.nworkers > 1 else ""
    pd.DataFrame(cav_rows).to_csv("%s/cavities_%s%s.csv" % (a.outdir, a.tag, sfx), index=False)
    if cell_frames:
        pd.concat(cell_frames, ignore_index=True).to_csv(
            "%s/cells_%s%s.csv" % (a.outdir, a.tag, sfx), index=False)
    naw = sum(1 for r in cav_rows if r["call"] == "airway")
    print("[%s w%d] DONE %.0f s | %d tiles | %d cavities | %d airways"
          % (a.tag, a.worker, time.time() - t0, done, len(cav_rows), naw), flush=True)


if __name__ == "__main__":
    main()
