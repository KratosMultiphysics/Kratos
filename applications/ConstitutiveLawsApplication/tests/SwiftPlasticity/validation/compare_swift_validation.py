"""Compare raw FE results against independent uniaxial bisection (no Kratos import)."""

from bisect import bisect_right
import csv
import hashlib
import importlib.util
import json
import math
import re
import struct

from validation_utils import DIRECTORY, read_reference_parameters, read_results, write_csv


def uniaxial_solution(strain, young_modulus, flow_stress):
    """Monotonic positive uniaxial stress: epsilon = sigma/E + p.

    Associated J2 flow gives d(epsilon_p_xx) = dp under uniaxial tension.
    Therefore F(p) = E*(epsilon-p) - sigma_y(p), with F strictly decreasing.
    In the plastic regime F(0)>0 and F(epsilon)<0 bracket the unique root.
    Bisection in total p is independent of the production incremental return
    mapping, plastic-multiplier normalization, and consistent FE tangent.
    """
    if young_modulus * strain <= flow_stress(0.0):
        return young_modulus * strain, 0.0
    lower, upper = 0.0, strain
    for _ in range(100):
        middle = 0.5 * (lower + upper)
        if middle == lower or middle == upper:
            break
        if young_modulus * (strain - middle) > flow_stress(middle):
            lower = middle
        else:
            upper = middle
    p = 0.5 * (lower + upper)
    sigma = young_modulus * (strain - p)
    # This checks the independent root's numerical accuracy, not FE agreement.
    if abs(sigma - flow_stress(p)) > 1e-9:
        raise RuntimeError("Analytical bisection did not satisfy its scalar equation")
    return sigma, p


def error_metrics(differences, relative_errors):
    meaningful = [value for value in relative_errors if value is not None]
    return {
        "max_absolute": max(abs(value) for value in differences),
        "rms_absolute": math.sqrt(math.fsum(value**2 for value in differences) / len(differences)),
        "max_relative": max(meaningful) if meaningful else None,
        "relative_error_rows": len(meaningful),
    }


def miso_diagnostic(parameters, ansys):
    """Quantify table interpolation separately; never use this for the Kratos material."""
    text = (DIRECTORY / "ansys/swift_validation_ansys.inp").read_text()
    table = [(float(p), float(sigma)) for p, sigma in re.findall(
        r"^TBPT,DEFI,([\d.eE+-]+),([\d.eE+-]+)\s*$", text, re.MULTILINE
    )]
    if len(table) != 100 or any(b[0] <= a[0] for a, b in zip(table, table[1:])):
        raise ValueError("Expected 100 strictly increasing APDL MISO plastic-strain points")
    abscissae = [point[0] for point in table]

    def tabulated_flow_stress(p):
        if not table[0][0] <= p <= table[-1][0]:
            raise ValueError("MISO diagnostic requires a root inside the supplied table")
        index = min(bisect_right(abscissae, p) - 1, len(table) - 2)
        p0, s0 = table[index]
        p1, s1 = table[index + 1]
        return s0 + (p - p0) * (s1 - s0) / (p1 - p0)

    differences = [row["sigma_xx"] - uniaxial_solution(
        row["eps_xx"], parameters["EMOD"], tabulated_flow_stress,
    )[0] for row in ansys]
    return {
        "points": len(table), "p_min": table[0][0], "p_max": table[-1][0],
        "max_table_ordinate_error_from_exact_swift": max(
            abs(sigma - parameters["KSWIFT"] * (parameters["EPS0"] + p)**parameters["NSWIFT"])
            for p, sigma in table
        ),
        "ansys_vs_tabulated_uniaxial_max_absolute_sigma": max(map(abs, differences)),
        "ansys_vs_tabulated_uniaxial_rms_absolute_sigma": math.sqrt(
            math.fsum(value**2 for value in differences) / len(differences)
        ),
        # This is a diagnostic of the supplied numbers, not a change to them or
        # a claim about ANSYS's internal arithmetic/storage implementation.
        "max_ansys_distance_to_float32": {
            key: max(abs(row[key] - struct.unpack("f", struct.pack("f", row[key]))[0])
                     for row in ansys) for key in ("sigma_xx", "seqv", "p_eq")
        },
    }


def make_plots(rows):
    if importlib.util.find_spec("matplotlib") is None:
        print("Matplotlib unavailable: plots skipped; no dependency installed.")
        return []
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = []
    for quantity, label, filename in (
        ("sigma", "Axial stress [MPa]", "swift_validation_stress_strain.png"),
        ("p", "Accumulated equivalent plastic strain p [-]",
         "swift_validation_plastic_strain.png"),
    ):
        figure, axes = plt.subplots(figsize=(7, 4.5))
        for source, legend, style in (
            ("analytical", "Analytical Swift", "-"),
            ("ansys", "ANSYS R19", "--"), ("kratos", "Kratos", ":"),
        ):
            axes.plot([row["eps_xx"] for row in rows],
                      [row[f"{quantity}_{source}"] for row in rows],
                      style, label=legend)
        axes.set(xlabel="Axial engineering strain [-]", ylabel=label,
                 title="One-element uniaxial Swift validation")
        axes.grid(True, alpha=0.3)
        axes.legend()
        figure.tight_layout()
        figure.savefig(DIRECTORY / filename, dpi=180)
        plt.close(figure)
        paths.append(filename)
    return paths


def main():
    parameters = read_reference_parameters()
    raw_paths = [DIRECTORY / name for name in (
        "ansys/swift_validation_ansys.inp", "ansys/swift_validation_ansys.csv",
        "kratos/swift_validation_kratos.csv",
    )]
    hashes = {str(path.relative_to(DIRECTORY)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in raw_paths}
    ansys = read_results(raw_paths[1], parameters)
    kratos = read_results(raw_paths[2], parameters)

    def swift_stress(p):
        return parameters["KSWIFT"] * (parameters["EPS0"] + p)**parameters["NSWIFT"]

    comparisons = ("ansys_analytical", "kratos_analytical", "kratos_ansys")
    # Relative errors use the second source as denominator. Values below these
    # physical thresholds are left blank in the CSV (null in JSON), not divided.
    thresholds = {"sigma": 1.0, "seqv": 1.0, "p": 1e-8}
    rows = []
    for ansys_row, kratos_row in zip(ansys, kratos):
        strain = ansys_row["eps_xx"]
        sigma, p = uniaxial_solution(strain, parameters["EMOD"], swift_stress)
        row = {"substep": int(ansys_row["substep"]), "eps_xx": strain}
        for quantity, source_column, analytical in (
            ("sigma", "sigma_xx", sigma), ("seqv", "seqv", abs(sigma)),
            ("p", "p_eq", p),
        ):
            values = {"analytical": analytical, "ansys": ansys_row[source_column],
                      "kratos": kratos_row[source_column]}
            row.update({f"{quantity}_{name}": value for name, value in values.items()})
            for comparison in comparisons:
                first, second = comparison.split("_")
                difference = values[first] - values[second]
                denominator = abs(values[second])
                suffix = f"{quantity}_{comparison}"
                row[f"difference_{suffix}"] = difference
                row[f"abs_error_{suffix}"] = abs(difference)
                row[f"relative_error_{suffix}"] = (
                    abs(difference) / denominator if denominator >= thresholds[quantity] else None
                )
        for source, solver_row in (("ansys", ansys_row), ("kratos", kratos_row)):
            for component in ("sigma_yy", "sigma_zz"):
                row[f"{component}_{source}"] = solver_row[component]
        rows.append(row)

    with (DIRECTORY / "kratos/swift_validation_diagnostics.csv").open(newline="") as stream:
        diagnostics = [{key: float(value) for key, value in row.items()}
                       for row in csv.DictReader(stream)]
    with (DIRECTORY / "kratos/swift_validation_integration_points.csv").open(newline="") as stream:
        points = [{key: float(value) for key, value in row.items()}
                  for row in csv.DictReader(stream)]
    if len(diagnostics) != 100 or len(points) != 800:
        raise ValueError("Expected 100 converged-step diagnostics and 800 integration-point rows")
    point_quantities = [key for key in points[0] if key not in ("substep", "integration_point")]
    maximum_spreads = dict.fromkeys(point_quantities, 0.0)
    for step in range(1, 101):
        step_points = [row for row in points if row["substep"] == step]
        if [row["integration_point"] for row in step_points] != list(range(1, 9)):
            raise ValueError(f"Missing or repeated integration points at step {step}")
        for key in point_quantities:
            maximum_spreads[key] = max(maximum_spreads[key],
                                      max(row[key] for row in step_points)
                                      - min(row[key] for row in step_points))

    summary = {
        "parameters_from_apdl": parameters, "raw_file_sha256": hashes,
        "relative_denominator_thresholds": thresholds,
        "initial_yield_stress": swift_stress(0.0),
        "initial_yield_strain": swift_stress(0.0) / parameters["EMOD"],
        "converged_increments": len(diagnostics), "integration_points": 8,
        "max_nonlinear_iterations": max(row["nonlinear_iterations"] for row in diagnostics),
        "max_final_residual_norm": max(row["residual_norm"] for row in diagnostics),
        "max_integration_point_spread": maximum_spreads,
        "errors": {
            quantity: {comparison: error_metrics(
                [row[f"difference_{quantity}_{comparison}"] for row in rows],
                [row[f"relative_error_{quantity}_{comparison}"] for row in rows],
            ) for comparison in comparisons} for quantity in thresholds
        },
        "max_transverse_stress": {
            source: {component: max(abs(row[component]) for row in source_rows)
                     for component in ("sigma_yy", "sigma_zz")}
            for source, source_rows in (("ansys", ansys), ("kratos", kratos))
        },
        "max_kratos_integration_point_transverse_stress": {
            component: max(abs(row[component]) for row in points)
            for component in ("sigma_yy", "sigma_zz")
        },
        "miso_diagnostic": miso_diagnostic(parameters, ansys),
        "final_step": {key: rows[-1][key] for key in (
            "eps_xx", "sigma_analytical", "sigma_ansys", "sigma_kratos",
            "seqv_analytical", "seqv_ansys", "seqv_kratos", "p_analytical", "p_ansys", "p_kratos",
        )},
    }
    write_csv(DIRECTORY / "swift_validation_comparison.csv", rows)
    summary["plots"] = make_plots(rows)
    for path in raw_paths:
        original_hash = hashes[str(path.relative_to(DIRECTORY))]
        if hashlib.sha256(path.read_bytes()).hexdigest() != original_hash:
            raise RuntimeError(f"Raw reference/result file changed during comparison: {path}")
    (DIRECTORY / "swift_validation_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
