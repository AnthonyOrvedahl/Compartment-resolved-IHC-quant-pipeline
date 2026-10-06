# Methods and parameters

Everything here is read from the code. Each value is a constant in the file named.

## 1. Data

Mouse lung sections, SARS-N immunohistochemistry (DAB, hematoxylin counterstain),
scanned on a Hamamatsu NanoZoomer C9600-12 (NDP.scan 3.1.9) at 20x, 0.4605 um/px;
ten slides scanned in one session. Slide #001 was used for development.

## 2. Stage 1: QuPath (`qupath/batch_classify_export.groovy`)

**Stain vectors** (optical density, background 255 255 255), estimated by Macenko
unmixing from a test region of slide #001 and applied to every slide:

| Stain | R | G | B |
|---|---|---|---|
| Hematoxylin | 0.71010 | 0.64228 | 0.28852 |
| DAB | 0.43668 | 0.58566 | 0.68287 |

**Cell detection** (`WatershedCellDetection`, hematoxylin OD image): pixel size 0.5 um,
background radius 8 um (by reconstruction), median radius 0, sigma 1.5 um, nucleus
area 10-400 um^2, threshold 0.1, max background 2.0, split by shape on, cell
expansion 7 um, smooth boundaries off, exclude DAB off. Detection runs inside the
existing annotation(s); a slide with two sections of the same block is annotated
around one section so tissue is not counted twice.

**Pixel classifier** (`Lung feature detection`): QuPath OpenCV random forest, input
resolution 0.4605 um/px, three channels (R, G, B), one feature per channel
(Gaussian, sigma 1.0), eleven classes: Airway lumen, Alveolar_Septa, Epithelial,
Epithelial wall, Indeterminate, Lymphocytic_infiltrate, Pleural_surface, Vessel
lumen, adipose, alveolar lumen, vessel wall. It sees colour at a single scale, so it
cannot separate airway from alveolar epithelium.

**Classification**: each detection takes the pixel classifier's class at its centroid.

## 3. Stage 2: geometry (`run_cohort.py`, `airway_acceptance.py`)

Per slide, tiles of 1400 um with a 400 um border read on every side, laid over the
bounding box of the exported cells.

1. Gaussian (sigma 1 px) on R, G, B at level 0; classifier applied every 3rd pixel
   (stride 3 reproduced the full-resolution prominence field with correlation 0.997
   and ran 9.2x faster); class map registered to the level-1 grid (0.921 um/px).
2. Tissue mask at level 1: hematoxylin OD + DAB OD > 0.15, deconvolved with the
   stain vectors above.
3. Class prominence: fraction of local tissue in a class group, smoothed at 10 um
   (adipose for the veto, 30 um), numerator and tissue mask smoothed identically
   and divided. Groups: epithelial = Epithelial + Epithelial wall; septal =
   Alveolar_Septa + alveolar lumen; lymphoid = Lymphocytic_infiltrate; adipose.
4. Candidate cavities: connected non-tissue regions of 200 um^2 up to 10% of the
   read window (tile plus border). A cavity is dropped if fewer than 40 tissue pixels lie in its ring.
5. Acceptance (constants in `PARAMS`), in order:

| Test | Rule | Value |
|---|---|---|
| adipose veto | median adipose prominence within 60 um > | 0.10 |
| lymphoid veto | median lymphoid prominence in the ring (2-20 um outside) > | 0.15 |
| airway | median(epithelial prominence - max(septal, lymphoid, adipose)) in the ring > | 0 |
| otherwise | rejected: septal competition | |

6. Each exported cell receives the distance from its centroid to the nearest
   accepted-airway cavity (`d_airway_um`).

Tiles are independent: a structure crossing a tile edge is scored once per side.

## 4. Stage 3: endpoints (`aggregate_cohort.py`)

| Rule | Value |
|---|---|
| DAB-positive | `Cell: DAB OD mean` >= 0.19161676646706588, one threshold for the whole cohort |
| Airway epithelium | classified Epithelial and `d_airway_um` <= 30 um |
| Alveolar+consolidated epithelium | all other Epithelial cells |
| Other reported compartments | Alveolar_Septa, Lymphocytic_infiltrate |
| Excluded before any denominator | Airway lumen, Vessel lumen, alveolar lumen, Indeterminate |
| Not reported (kept in the infiltrate denominator) | Epithelial wall, vessel wall, Pleural_surface, adipose |
| Unit of analysis | the mouse; no minimum-cell floor, precision shown as a Wilson 95% interval |
| Group comparison | two-sided Mann-Whitney per compartment and day, and pooled over days |
| Compartment comparison | paired Wilcoxon on all mice with every compartment, Holm-adjusted; Friedman test across the four |
| Attainable p | `min_possible_p` is reported; when above 0.05 a result cannot reach significance at that group size |

All p-values are unadjusted except `p_holm`. With 2-3 mice per group and day, the
smallest attainable two-sided p is 0.10 (3 v 3) and 0.33 (2 v 2).
Infiltrate fraction = Lymphocytic_infiltrate cells / all analysed cells per mouse.

## 5. Basis of the constants

| Constant | Basis |
|---|---|
| Stain vectors | Macenko unmixing, one region of slide #001 |
| DAB threshold | the investigator's QuPath object classifier (`config/DAB_pos_cells.json`); no negative-control slide was available, so it was set from image appearance |
| Acceptance constants | one region of slide #001: margin sign from 5 confirmed airways (+0.33..+0.69) against 3 alveolar cavities (-0.21..-0.33); lymphoid veto fitted to 7 investigator calls; adipose veto separates airway calls (max 0.050) from fat vacuoles (0.143-0.206) |
| 30 um airway band | epithelial thickness measured on this material, 18-27 um |
| Tile 1400 / border 400 / stride 3 | tile size bounds memory (one read is 2200 x 2200 um); the 400 um border lets structures crossing a tile edge be seen; stride from the benchmark above |

## 6. Limitations

* Airway acceptance was set and checked on one region of one slide. It has not been
  validated on the other nine.
* Merged voids: a bronchus lumen joined to the alveolar airspace through gaps in the
  tissue mask forms one large cavity, whose ring samples alveolar septum, so it is
  rejected. Airways of this kind are missed.
* The pixel classifier reads R, G and B and so sees DAB. Airway acceptance has not
  been tested for dependence on staining intensity.
* The QuPath class Epithelial includes epithelium in consolidated tissue; the airway
  split rests only on distance to an accepted airway.
* DAB OD depends on the stain vectors set on each image. Run `check_stains.py` on
  every slide.
* The classifier's training annotations are not distributed.
