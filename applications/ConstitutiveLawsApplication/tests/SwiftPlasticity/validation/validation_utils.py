"""Shared material input, load grid and output helpers for the FE benchmark."""

import argparse
import csv
import json
import math
from pathlib import Path
import tempfile


DIRECTORY = Path(__file__).resolve().parent
MATERIAL_PATH = DIRECTORY / "kratos/Materials.json"
LENGTH = 1.0  # mm
MAX_DISPLACEMENT = 0.02  # mm
NUMBER_OF_STEPS = 100
COLUMNS = (
    "substep", "eps_xx", "sigma_xx", "sigma_yy", "sigma_zz",
    "seqv", "p_eq", "ux",
)


def parse_output_directory(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(tempfile.gettempdir()) / "kratos_swift_validation",
        help="Result directory (default: %(default)s); reused on reruns.",
    )
    return parser.parse_args().output_dir.expanduser().resolve()


def read_material_parameters():
    properties = json.loads(MATERIAL_PATH.read_text())["properties"][0]
    material = properties["Material"]
    return material["Variables"]


def read_results(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames) != COLUMNS:
            raise ValueError(f"Unexpected columns in {path}")
        rows = [{key: float(value) for key, value in row.items()}
                for row in reader]
    count = NUMBER_OF_STEPS
    if len(rows) != count:
        raise ValueError(
            f"Expected {count} result rows in {path}, got {len(rows)}",
        )
    for step, row in enumerate(rows, 1):
        if not all(math.isfinite(value) for value in row.values()):
            raise ValueError(f"Non-finite result at step {step} in {path}")
        expected_ux = MAX_DISPLACEMENT * step / count
        if (row["substep"] != step
                or abs(row["ux"] - expected_ux) > 1e-14
                or abs(row["eps_xx"] - expected_ux / LENGTH) > 1e-14):
            raise ValueError(f"Unexpected load grid at step {step} in {path}")
    return rows


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
