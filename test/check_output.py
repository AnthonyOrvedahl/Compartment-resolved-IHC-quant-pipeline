#!/usr/bin/env python
"""Compare run_cohort.py / aggregate_cohort.py test output with test/expected/.

Usage: python test/check_output.py <output-dir> [--stats]

Exact: same cells, same candidate cavities, same airway/rejected calls.
Tolerant: floating-point columns to 1e-3 (distances in um, margins), because
OpenCV, SciPy and NumPy builds can differ in the last digits.
Exit status 0 = pass.
"""
import os
import sys

import numpy as np
import pandas as pd

here = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1]
ok = True


def compare(name, key=None, tol=1e-3):
    global ok
    a = pd.read_csv(os.path.join(out, name), dtype={"slide": str})
    b = pd.read_csv(os.path.join(here, "expected", name), dtype={"slide": str})
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        print("FAIL %s: shape/columns differ %s vs %s" % (name, a.shape, b.shape))
        ok = False
        return
    if key:
        a, b = a.sort_values(key).reset_index(drop=True), b.sort_values(key).reset_index(drop=True)
    bad = []
    for c in b.columns:
        if pd.api.types.is_float_dtype(b[c]):
            if not np.allclose(a[c], b[c], atol=tol, rtol=0, equal_nan=True):
                bad.append(c)
        elif not (a[c].astype(str) == b[c].astype(str)).all():
            bad.append(c)
    print("%s %s%s" % ("PASS" if not bad else "FAIL", name, "" if not bad else "  columns differ: %s" % bad))
    ok = ok and not bad


if "--stats" in sys.argv:
    compare("group_contrasts.csv", tol=1e-4)
    compare("compartment_contrasts.csv", tol=1e-4)
else:
    compare("cells_test.csv", key="Object ID")
    compare("cavities_test.csv", key="cavity")
sys.exit(0 if ok else 1)
