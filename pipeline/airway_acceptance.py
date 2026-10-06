"""Airway-lumen acceptance from pixel-classifier class prominences.

Input: a hard-class map from the pixel classifier, a tissue mask, the pixel size.
All lengths are micrometres; every constant is in PARAMS.

Steps
  1. Class prominence = fraction of local tissue belonging to a class, smoothed
     at prom_sigma_um (adipose for the veto: adi_sigma_um). The numerator and the
     tissue mask are smoothed identically and divided; without the division the
     airspace inside a lumen dilutes the epithelial signal at exactly the
     structures of interest.
  2. Candidate cavities = connected non-tissue regions larger than cav_min_um2 and
     smaller than cav_max_frac of the field.
  3. For each cavity, in the ring ring_inner_um..ring_outer_um outside it:
         margin = median(epithelial prominence - max(septal, lymphoid, adipose))
  4. Call, in this order:
         adipose veto    median adipose prominence within adi_ring_um > adi_veto
         lymphoid veto   median lymphoid prominence in the ring     > lym_veto
         airway          margin > margin_min
         otherwise       rejected: septal competition

Provenance of the constants (development slide #001, one region; see METHODS.md)
  margin_min = 0   five investigator-confirmed airways scored +0.33..+0.69 and three
                   alveolar cavities -0.21..-0.33 (gap 0.54, no overlap), so the sign
                   of zero is the threshold. Not tuned.
  lym_veto  = 0.15 fitted to seven investigator calls in one region.
  adi_veto  = 0.10 airway calls reach at most 0.050; confirmed fat vacuoles
                   0.143-0.206.
  cav_min_um2 = 200 compressed distal airways have small lumens; at 800 a confirmed
                   airway with slit voids of 329 and 309 um^2 was never segmented.

Known limitations
  * Merged voids: a bronchus lumen connected to the alveolar airspace through gaps
    in the tissue mask is one large cavity whose ring samples mostly alveolar
    septum, so it is rejected. Airways of this kind are missed.
  * The class maps come from a classifier that reads R, G and B and therefore sees
    DAB. Airway acceptance has not been tested for dependence on staining intensity.
  * lym_veto and adi_veto were set on a single region.
"""

import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter
from skimage.morphology import remove_small_objects

CLASS_GROUPS = {
    "EPI": ["Epithelial", "Epithelial wall"],
    "SEPT": ["Alveolar_Septa", "alveolar lumen"],
    "LYM": ["Lymphocytic_infiltrate"],
    "ADI": ["adipose"],
}

PARAMS = dict(
    prom_sigma_um=10.0,      # tissue scale for class prominence
    adi_sigma_um=30.0,       # wider context for the adipose veto
    ring_inner_um=2.0,
    ring_outer_um=20.0,
    adi_ring_um=60.0,
    cav_min_um2=200.0,       # NOT 800 - see docstring
    cav_max_frac=0.10,
    lym_veto=0.15,           # provisional, one ROI, 7 labels
    adi_veto=0.10,
    margin_min=0.0,          # sign of zero; no tuning
)


def prominence(indicator, tissue, mpp, sigma_um):
    """Class fraction of LOCAL TISSUE at tissue scale.

    Dividing by the identically smoothed tissue mask is load-bearing: without
    it the airspace inside a lumen dilutes the numerator and the epithelial
    signal collapses at exactly the structures of interest.
    """
    num = gaussian_filter((indicator & tissue).astype(np.float32), sigma_um / mpp)
    den = gaussian_filter(tissue.astype(np.float32), sigma_um / mpp) + 1e-6
    return num / den


def class_prominences(pred, class_names, tissue, mpp, p=None):
    p = dict(PARAMS, **(p or {}))
    out = {}
    for key, names in CLASS_GROUPS.items():
        idx = [class_names.index(n) for n in names if n in class_names]
        ind = np.isin(pred, idx) if idx else np.zeros(pred.shape, bool)
        sig = p["adi_sigma_um"] if key == "ADI" else p["prom_sigma_um"]
        out[key] = prominence(ind, tissue, mpp, sig)
    out["ADI_ring"] = prominence(
        np.isin(pred, [class_names.index("adipose")]) if "adipose" in class_names
        else np.zeros(pred.shape, bool), tissue, mpp, p["prom_sigma_um"])
    return out


def enumerate_cavities(tissue, mpp, p=None):
    p = dict(PARAMS, **(p or {}))
    void = remove_small_objects(~tissue, max(int(p["cav_min_um2"] / mpp ** 2), 1))
    lab, _ = ndi.label(void)
    size = np.bincount(lab.ravel()) * mpp ** 2
    size[0] = 0
    field = tissue.size * mpp ** 2
    ok = (size > p["cav_min_um2"]) & (size < p["cav_max_frac"] * field)
    ok[0] = False
    return lab, np.where(ok)[0], size


def score_cavities(lab, ids, sizes, proms, tissue, mpp, p=None):
    """One row per candidate cavity with the margin, both vetoes and the call."""
    p = dict(PARAMS, **(p or {}))
    rival = np.maximum.reduce([proms["SEPT"], proms["LYM"], proms["ADI_ring"]])
    rows = []
    for cav in ids:
        m = lab == cav
        dist = ndi.distance_transform_edt(~m) * mpp
        ring = tissue & (dist > p["ring_inner_um"]) & (dist < p["ring_outer_um"])
        wide = tissue & (dist < p["adi_ring_um"])
        if ring.sum() < 40:
            continue
        margin = float(np.median(proms["EPI"][ring] - rival[ring]))
        lym = float(np.median(proms["LYM"][ring]))
        adi = float(np.median(proms["ADI"][wide])) if wide.any() else 0.0
        if adi > p["adi_veto"]:
            call = "rejected: adipose veto"
        elif lym > p["lym_veto"]:
            call = "rejected: lymphoid veto"
        elif margin > p["margin_min"]:
            call = "airway"
        else:
            call = "rejected: septal competition"
        rows.append(dict(cavity=int(cav), area_um2=float(sizes[cav]),
                         equiv_diam_um=float(2 * np.sqrt(sizes[cav] / np.pi)),
                         margin=margin, lym_ring=lym, adi_wide=adi, call=call))
    return rows
