# SARS-N IHC: compartment-resolved quantification in mouse lung

Measures SARS-CoV-2 nucleocapsid (DAB) positivity **per tissue compartment** in
whole lung sections: airway epithelium, alveolar and consolidated epithelium,
alveolar septa, lymphocytic infiltrate. Mice are Cas9-EGFP transgenic, transduced
with an AdV carrying hACE2-mKate2 and an sgRNA against EGFP (control) or *Tmprss2*,
infected with SARS-CoV-2 WA1 and harvested at 3 or 7 days post-infection. The
unit of analysis is the mouse.

> **Status.** Cohort result tables and figures are not included. Two problems found
> while preparing this package change them, and the cohort has to be re-run first:
> 1. `run_cohort.py` registered the classifier's class map to the tissue mask at the
>    wrong scale. Fixed here; the earlier code is not reproduced by this version.
> 2. Slide #002's exported DAB values were computed with QuPath's default stain
>    vectors; the other nine slides used the cohort vectors. Re-detect #002 with the
>    cohort vectors (`pipeline/check_stains.py` reports the mismatch), then re-run.
>
> Delete this block once the cohort has been re-run and `results/` added.

## How it works

| Stage | Script | Output |
|---|---|---|
| 1. Detect cells, classify each by tissue type, export | `qupath/batch_classify_export.groovy` (QuPath 0.7.0) | `exports/<image>.Classifications.txt` |
| QC. Stain vectors match the cohort vectors? | `pipeline/check_stains.py` | pass/fail per slide |
| 2. Find airways, give each cell its distance to the nearest one | `pipeline/run_cohort.py`, `pipeline/airway_acceptance.py` | `cells_<tag>.csv`, `cavities_<tag>.csv` |
| 3. Per-mouse endpoints and statistics | `pipeline/aggregate_cohort.py` | four tables + `summary.json` |

QuPath decides what each cell *is*; Python decides where it sits. Airway versus
alveolar epithelium are the same colour, so that distinction comes from geometry.
Rules, every constant and its basis are in [METHODS.md](METHODS.md).

All Python is deterministic: fixed rules and one fixed trained forest, no random
number generation, no model calls. Two runs on the same input give byte-identical output (tested).

## Requirements

Python 3.13 (tested). QuPath 0.7.0 is needed for stage 1 only.

    pip install -r requirements.txt

That is all: `openslide-bin` bundles the OpenSlide library, so nothing is installed
system-wide. This was tested in a fresh environment on macOS, where the package
tests below pass and reproduce the expected output byte for byte. If pip has no
wheel for your platform, install OpenSlide with your system package manager
(`brew install openslide`, `apt install libopenslide0`, or conda-forge) and delete
the `openslide-bin` line. Each `run_cohort.py` process needs about 2-3 GB of RAM.

## Run

    # 1. QuPath: open the project, Automate > Run for project > batch_classify_export.groovy
    # 2. QC (exit status 1 on mismatch)
    python pipeline/check_stains.py --slide S.ndpi --meas exports/S.Classifications.txt
    # 3. per slide; shard across workers if wanted
    for w in 0 1 2 3 4 5 6 7; do
      python pipeline/run_cohort.py --slide S.ndpi --meas exports/S.Classifications.txt \
          --classifier "Lung feature prediction.json" --outdir out --tag 001 \
          --worker $w --nworkers 8 &
    done; wait
    # 4. all slides listed in the manifest (columns: slide, sample, group, day)
    python pipeline/aggregate_cohort.py --results-dir out --manifest config/manifest.csv --out-dir results

The manifest in `config/manifest_example.csv` assumes slide order follows the row
order of the sample key; confirm that mapping before use.

## Test (about 3 minutes in total)

`test/` holds a 1.39 mm square cut from slide #001 and the matching QuPath export.

    # with Lung feature prediction.json (see below) in the package root
    python pipeline/run_cohort.py --slide test/crop_slide001.tif \
        --meas test/crop_slide001_cells.tsv --classifier "Lung feature prediction.json" \
        --outdir out_test --tag test
    python test/check_output.py out_test                       # expect PASS PASS
    python pipeline/aggregate_cohort.py --per-mouse test/stats_fixture_per_mouse.csv --out-dir out_stats
    python test/check_output.py out_stats --stats              # expect PASS PASS

The first test needs the classifier (below), runs one tile (14 airways expected)
and takes about 2.5 minutes on a 10-core laptop. The second needs nothing external and runs in under a second. Its input is a
fixture for the statistics code only; its numbers come from the pre-fix run and are not results.
`check_output.py` compares exactly for calls and counts and to 1e-3 for floats.

| File | Size | Purpose |
|---|---|---|
| `test/crop_slide001.tif` | 8.6 MB | tiled TIFF, 0.4605 um/px |
| `test/crop_slide001_cells.tsv` | 1.9 MB | 18,151 cells, six columns, coordinates shifted to the crop |
| `test/expected/*` | 1.5 MB | outputs of the tests above |

## Not in this repository

| Item | Size | Why / where |
|---|---|---|
| `Lung feature prediction.json`, sha256 `e022338361fe9cfda3028703180508f0af8b9b1c95b289dfab241c430ae978a8` | 374 MB | exceeds git limits; deposit externally. Used by `run_cohort.py` |
| `Lung feature detection.json`, sha256 `944674ccfdb332f477cf301c6c5079acc8aecf1ead9c13b17c5447cd54441ae6` | 374 MB | same trained forest, classification output; load into the QuPath project for stage 1 |
| 10 NDPI slides | 1.6 GB | raw data |
| 10 QuPath exports | about 0.8 GB | regenerate with stage 1 |

The images used to train the pixel classifier are not included, so it can be used
but not retrained from this repository.

## Layout

    pipeline/   run_cohort.py  airway_acceptance.py  aggregate_cohort.py  check_stains.py
    qupath/     batch_classify_export.groovy
    config/     manifest_example.csv  DAB_pos_cells.json (source of the positivity threshold)
    test/       crop, export, expected outputs, check_output.py

The Groovy script was not run while preparing this package (QuPath was unavailable);
its only change from the version used is the output directory, now `<project>/exports`.

Add a LICENSE before publishing.
