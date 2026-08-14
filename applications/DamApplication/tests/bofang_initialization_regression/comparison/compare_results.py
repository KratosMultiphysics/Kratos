"""Compare the Bofang initialization-regression variants.

Variants (results stored under ``results/<name>``):
    A : current_master
    B : current_master_legacy_bofang
    C : pre_13472

The most important comparison is A vs B.
"""

import argparse
import os

import pandas as pd


VARIANT_DIRS = {
    "A": "current_master",
    "B": "current_master_legacy_bofang",
    "B2": "current_master_legacy_bofang_faithful",
    "C": "pre_13472",
}


def load_results(base_dir, variant):
    results = {}
    for fname in ["lifecycle_temperatures.csv", "results.csv", "stress_results.csv"]:
        path = os.path.join(base_dir, VARIANT_DIRS[variant], fname)
        results[fname] = pd.read_csv(path)
    return results


def report_max(name, df_a, df_b, numeric_cols, key_cols, rel_cols=None):
    merged = df_a.merge(df_b, on=key_cols, suffixes=("_A", "_B"), how="outer")
    print("=" * 78)
    print(name)
    print("=" * 78)
    for col in numeric_cols:
        ca = col + "_A"
        cb = col + "_B"
        diff = (merged[ca] - merged[cb]).abs()
        idx = diff.idxmax()
        val = diff[idx]
        if rel_cols and col in rel_cols:
            denom = merged[ca].abs().max()
            rel = val / denom if denom > 0 else float("nan")
        else:
            rel = float("nan")
        row = merged.loc[idx, key_cols].to_dict()
        print("  %-22s max_abs_diff = %.3e   (rel to A max %.3e)  at %s"
              % (col, val, rel, row))


def compare_pair(base_dir, name_a, name_b):
    ra = load_results(base_dir, name_a)
    rb = load_results(base_dir, name_b)

    # lifecycle temperatures
    lifecycle_keys = ["lifecycle_stage", "node_id"]
    lc_cols = ["temperature"]
    report_max(
        "Lifecycle TEMPERATURE  %s vs %s" % (name_a, name_b),
        ra["lifecycle_temperatures.csv"], rb["lifecycle_temperatures.csv"],
        lc_cols, lifecycle_keys, rel_cols=lc_cols,
    )

    # step results
    results_keys = ["time", "node_id"]
    res_cols = ["temperature", "displacement_x", "displacement_y", "displacement_norm"]
    report_max(
        "Step results  %s vs %s" % (name_a, name_b),
        ra["results.csv"], rb["results.csv"], res_cols, results_keys, rel_cols=res_cols,
    )

    # stresses
    stress_keys = ["time", "node_id"]
    stress_cols = ["stress_xx", "stress_yy", "stress_xy"]
    report_max(
        "Stress results  %s vs %s" % (name_a, name_b),
        ra["stress_results.csv"], rb["stress_results.csv"],
        stress_cols, stress_keys, rel_cols=stress_cols,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", default="results")
    parser.add_argument("--pairs", default="A-B,B-C,A-C",
                        help="comma-separated variant pairs to compare (default: A-B,B-C,A-C)")
    args = parser.parse_args()

    for fname, variant in VARIANT_DIRS.items():
        for f in ["lifecycle_temperatures.csv", "results.csv", "stress_results.csv"]:
            path = os.path.join(args.base_dir, variant, f)
            if not os.path.exists(path):
                print("NOTE: missing expected result file: %s (variant skipped)" % path)
        if os.path.exists(os.path.join(args.base_dir, variant, "lifecycle_temperatures.csv")):
            print("Loaded variant %s (%s)" % (fname, variant))

    for pair in args.pairs.split(","):
        name_a, name_b = pair.split("-")
        da = os.path.join(args.base_dir, VARIANT_DIRS[name_a], "lifecycle_temperatures.csv")
        db = os.path.join(args.base_dir, VARIANT_DIRS[name_b], "lifecycle_temperatures.csv")
        if not (os.path.exists(da) and os.path.exists(db)):
            print()
            print("### SKIPPING %s vs %s (missing results)" % (name_a, name_b))
            continue
        print()
        print("### COMPARISON %s vs %s" % (name_a, name_b))
        compare_pair(args.base_dir, name_a, name_b)


if __name__ == "__main__":
    main()
