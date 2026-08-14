"""Analytical vs Kratos-assigned Bofang temperatures for every variant.

Compares the analytical Bofang expression (implemented in
comparison/analytical_bofang.py) against the temperature that Kratos assigned
to the reservoir nodes right after the full analysis initialization
(lifecycle stage "after_execute_before_solution_loop").
"""

import argparse
import os

from analytical_bofang import BofangAnalytical

VARIANT_DIRS = {
    "A": "current_master",
    "B": "current_master_legacy_bofang",
    "B2": "current_master_legacy_bofang_faithful",
    "C": "pre_13472",
}

FACE_NODE_Y = {1: 20.0, 4: 16.0, 7: 12.0, 10: 8.0, 13: 4.0, 16: 0.0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", default="../results")
    args = parser.parse_args()

    analytical = BofangAnalytical()
    stage = "after_execute_before_solution_loop"

    assigned = {}
    for variant, dname in VARIANT_DIRS.items():
        path = os.path.join(args.base_dir, dname, "lifecycle_temperatures.csv")
        values = {}
        if os.path.exists(path):
            with open(path) as f:
                header = f.readline().strip().split(",")
                for line in f:
                    rec = dict(zip(header, line.strip().split(",")))
                    if rec["lifecycle_stage"] == stage:
                        values[int(rec["node_id"])] = float(rec["temperature"])
        assigned[variant] = values

    print("node_id, elevation, depth, analytical, " +
          ", ".join("assigned_%s, abs_err_%s" % (v, v) for v in VARIANT_DIRS))
    for node_id in sorted(FACE_NODE_Y):
        y = FACE_NODE_Y[node_id]
        depth = analytical.water_level - y
        expected = analytical.temperature({"id": node_id, "x": 0.0, "y": y})
        row = [node_id, y, depth if depth >= 0 else "n/a",
               "%.10f" % expected if expected is not None else "n/a"]
        for variant in VARIANT_DIRS:
            val = assigned[variant].get(node_id)
            if expected is None:
                row += ["0.0" if val == 0.0 else "%.6e" % (val or 0.0), "-"]
            elif val is not None:
                row += ["%.10f" % val, "%.3e" % abs(val - expected)]
            else:
                row += ["n/a", "n/a"]
        print(", ".join(str(x) for x in row))


if __name__ == "__main__":
    main()
