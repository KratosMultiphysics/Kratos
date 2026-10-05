"""Read-only helpers for this specific APDL reference, not an APDL interpreter."""

import csv
import math
from pathlib import Path
import re


DIRECTORY = Path(__file__).resolve().parent
COLUMNS = (
    "substep", "eps_xx", "sigma_xx", "sigma_yy", "sigma_zz", "seqv", "p_eq", "ux"
)


def read_reference_parameters():
    text = (DIRECTORY / "ansys/swift_validation_ansys.inp").read_text()
    names = ("EMOD", "NU", "KSWIFT", "EPS0", "NSWIFT", "LENGTH", "UMAX", "NSUB")
    parameters = {}
    for name in names:
        match = re.search(rf"^{name}\s*=\s*([\d.eE+-]+)\s*$", text, re.MULTILINE)
        if match is None:
            raise ValueError(f"Missing numeric APDL assignment: {name}")
        parameters[name] = float(match[1])
    if parameters["NSUB"] != 100 or parameters["UMAX"] / parameters["LENGTH"] != 0.02:
        raise ValueError("This benchmark requires 100 increments ending at strain 0.02")
    return parameters


def read_results(path, parameters):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames) != COLUMNS:
            raise ValueError(f"Unexpected columns in {path}")
        rows = [{key: float(value) for key, value in row.items()} for row in reader]
    count = int(parameters["NSUB"])
    if len(rows) != count:
        raise ValueError(f"Expected {count} result rows in {path}, got {len(rows)}")
    for step, row in enumerate(rows, 1):
        if not all(math.isfinite(value) for value in row.values()):
            raise ValueError(f"Non-finite result at step {step} in {path}")
        expected_ux = parameters["UMAX"] * step / count
        if (row["substep"] != step
                or abs(row["ux"] - expected_ux) > 1e-14
                or abs(row["eps_xx"] - expected_ux / parameters["LENGTH"]) > 1e-14):
            raise ValueError(f"Unexpected load grid at step {step} in {path}")
    return rows


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
