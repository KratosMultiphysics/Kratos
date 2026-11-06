"""Compare FE results with independent bisection, without importing Kratos."""

import csv
import hashlib
import importlib.util
import json
import math

from validation_utils import (
    NUMBER_OF_STEPS, parse_output_directory, read_material_parameters,
    read_results, write_csv,
)


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
    # Check the independent root's numerical accuracy, not FE agreement.
    if abs(sigma - flow_stress(p)) > 1e-9:
        raise RuntimeError(
            "Analytical bisection did not satisfy its scalar equation",
        )
    return sigma, p


def error_metrics(differences, relative_errors):
    meaningful = [value for value in relative_errors if value is not None]
    return {
        "max_absolute": max(abs(value) for value in differences),
        "rms_absolute": math.sqrt(
            math.fsum(value**2 for value in differences) / len(differences)
        ),
        "max_relative": max(meaningful) if meaningful else None,
        "relative_error_rows": len(meaningful),
    }


def make_plots(rows, output):
    if importlib.util.find_spec("matplotlib") is None:
        print("Matplotlib unavailable: skipping plots.")
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
            ("kratos", "Kratos FE", ":"),
        ):
            axes.plot([row["eps_xx"] for row in rows],
                      [row[f"{quantity}_{source}"] for row in rows],
                      style, label=legend)
        axes.set(xlabel="Axial engineering strain [-]", ylabel=label,
                 title="One-element uniaxial Swift validation")
        axes.grid(True, alpha=0.3)
        axes.legend()
        figure.tight_layout()
        figure.savefig(output / filename, dpi=180)
        plt.close(figure)
        paths.append(filename)
    return paths


def compare_results(output):
    """Read the raw FE outputs and write analytical comparisons beside them."""
    parameters = read_material_parameters()
    raw_path = output / "swift_validation_kratos.csv"
    original_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    kratos = read_results(raw_path)
    young_modulus = parameters["YOUNG_MODULUS"]

    def swift_stress(p):
        return parameters["SWIFT_COEFFICIENT"] * (
            parameters["SWIFT_INITIAL_STRAIN"] + p
        )**parameters["SWIFT_HARDENING_EXPONENT"]

    # Use analytical magnitudes as denominators. Below these physical
    # thresholds leave CSV cells blank (JSON null); errors are fractions.
    thresholds = {"sigma": 1.0, "seqv": 1.0, "p": 1e-8}
    rows = []
    for kratos_row in kratos:
        strain = kratos_row["eps_xx"]
        sigma, p = uniaxial_solution(strain, young_modulus, swift_stress)
        row = {"substep": int(kratos_row["substep"]), "eps_xx": strain}
        for quantity, source_column, analytical in (
            ("sigma", "sigma_xx", sigma), ("seqv", "seqv", abs(sigma)),
            ("p", "p_eq", p),
        ):
            value = kratos_row[source_column]
            difference = value - analytical
            denominator = abs(analytical)
            row[f"{quantity}_analytical"] = analytical
            row[f"{quantity}_kratos"] = value
            row[f"difference_{quantity}"] = difference
            row[f"abs_error_{quantity}"] = abs(difference)
            row[f"relative_error_{quantity}"] = (
                abs(difference) / denominator
                if denominator >= thresholds[quantity] else None
            )
        for component in ("sigma_yy", "sigma_zz"):
            row[component] = kratos_row[component]
        rows.append(row)

    diagnostic_path = output / "swift_validation_diagnostics.csv"
    with diagnostic_path.open(newline="") as stream:
        diagnostics = [{key: float(value) for key, value in row.items()}
                       for row in csv.DictReader(stream)]
    point_path = output / "swift_validation_integration_points.csv"
    with point_path.open(newline="") as stream:
        points = [{key: float(value) for key, value in row.items()}
                  for row in csv.DictReader(stream)]
    if (len(diagnostics) != NUMBER_OF_STEPS
            or len(points) != 8 * NUMBER_OF_STEPS):
        raise ValueError(
            "Expected 100 converged-step diagnostics and 800 point rows",
        )
    point_quantities = [key for key in points[0]
                        if key not in ("substep", "integration_point")]
    maximum_spreads = dict.fromkeys(point_quantities, 0.0)
    for step in range(1, NUMBER_OF_STEPS + 1):
        step_points = [row for row in points if row["substep"] == step]
        point_ids = [row["integration_point"] for row in step_points]
        if point_ids != list(range(1, 9)):
            raise ValueError(f"Missing or repeated points at step {step}")
        for key in point_quantities:
            spread = (max(row[key] for row in step_points)
                      - min(row[key] for row in step_points))
            maximum_spreads[key] = max(maximum_spreads[key], spread)

    summary = {
        "material_parameters": parameters,
        "raw_result_sha256": original_hash,
        "relative_denominator_thresholds": thresholds,
        "initial_yield_stress": swift_stress(0.0),
        "initial_yield_strain": swift_stress(0.0) / young_modulus,
        "converged_increments": len(diagnostics), "integration_points": 8,
        "max_nonlinear_iterations": max(
            row["nonlinear_iterations"] for row in diagnostics
        ),
        "max_final_residual_norm": max(
            row["residual_norm"] for row in diagnostics
        ),
        "max_integration_point_spread": maximum_spreads,
        "errors": {
            quantity: error_metrics(
                [row[f"difference_{quantity}"] for row in rows],
                [row[f"relative_error_{quantity}"] for row in rows],
            ) for quantity in thresholds
        },
        "max_transverse_stress": {
            component: max(abs(row[component]) for row in kratos)
            for component in ("sigma_yy", "sigma_zz")
        },
        "max_integration_point_transverse_stress": {
            component: max(abs(row[component]) for row in points)
            for component in ("sigma_yy", "sigma_zz")
        },
        "final_step": {key: rows[-1][key] for key in (
            "eps_xx", "sigma_analytical", "sigma_kratos", "seqv_analytical",
            "seqv_kratos", "p_analytical", "p_kratos",
        )},
    }
    write_csv(output / "swift_validation_comparison.csv", rows)
    summary["plots"] = make_plots(rows, output)
    if hashlib.sha256(raw_path.read_bytes()).hexdigest() != original_hash:
        raise RuntimeError(f"Raw FE result changed: {raw_path}")
    (output / "swift_validation_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    compare_results(parse_output_directory(__doc__))
