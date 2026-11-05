"""Run the structured steel-plate analyses against the local responses.py."""

import argparse
import csv
import importlib.util
import json
from math import isclose, sqrt
from pathlib import Path

import KratosMultiphysics as KM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)


ROOT = Path(__file__).resolve().parent
CASES = (
    "plane_stress_q4", "plane_stress_q8",
    "plane_strain_q4", "plane_strain_q8",
    "solid_h8", "solid_h20",
)


def create_analysis(case_directory):
    """Resolve input paths so Main.py can be launched from any directory."""
    parameters = KM.Parameters((case_directory / "ProjectParameters.json").read_text())
    settings = parameters["solver_settings"]
    settings["model_import_settings"]["input_filename"].SetString(
        str(case_directory / "plate")
    )
    settings["material_import_settings"]["materials_filename"].SetString(
        str(case_directory / "Materials.json")
    )
    model = KM.Model()
    return model, StructuralMechanicsAnalysis(model, parameters)


def load_response_extractor():
    """Load the working-tree helper, even when Kratos uses an installed build."""
    source = ROOT.parents[1] / "python_scripts" / "helper" / "responses.py"
    spec = importlib.util.spec_from_file_location("local_structural_responses", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    print(f"Response source: {source}")
    return module.ResponseExtractor()


def write_results(case_directory, model_part, values, mode):
    results_directory = case_directory / "results"
    results_directory.mkdir(exist_ok=True)
    downstream_values = []
    with (results_directory / f"von_mises_{mode}.csv").open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["element_id", "center_x_mm", "center_y_mm", "center_z_mm",
                         "gauss_points", "von_mises_mpa"])
        for element, value in zip(model_part.Elements, values):
            center = element.GetGeometry().Center()
            writer.writerow([element.Id, center.X, center.Y, center.Z,
                             len(element.GetIntegrationPoints()), float(value)])
            if center.X >= 75.0:
                downstream_values.append(float(value))
    with (results_directory / "nodal_results.csv").open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["node_id", "x_mm", "y_mm", "z_mm", "ux_mm", "uy_mm", "uz_mm",
                         "reaction_x_n", "reaction_y_n", "reaction_z_n"])
        for node in model_part.Nodes:
            writer.writerow([node.Id, node.X0, node.Y0, node.Z0,
                             *node.GetSolutionStepValue(KM.DISPLACEMENT),
                             *node.GetSolutionStepValue(KM.REACTION)])
    reaction_x = sum(node.GetSolutionStepValue(KM.REACTION_X)
                     for node in model_part.GetSubModelPart("Fixed").Nodes)
    reference = 100.0 * sqrt(1.0 - 0.3 + 0.3**2) if case_directory.name.startswith("plane_strain") else 100.0
    summary = {
        "case": case_directory.name,
        "response_mode": mode,
        "nodes": model_part.NumberOfNodes(),
        "elements": model_part.NumberOfElements(),
        "von_mises_min_mpa": float(min(values)),
        "von_mises_max_mpa": float(max(values)),
        "downstream_mean_von_mises_mpa": sum(downstream_values) / len(downstream_values),
        "ideal_far_field_von_mises_mpa": reference,
        "loaded_edge_mean_ux_mm": sum(
            node.GetSolutionStepValue(KM.DISPLACEMENT_X)
            for node in model_part.GetSubModelPart("Loaded").Nodes
        ) / model_part.GetSubModelPart("Loaded").NumberOfNodes(),
        "applied_force_x_n": 4000.0,
        "support_reaction_x_n": reaction_x,
        "force_balance_error_n": reaction_x + 4000.0,
    }
    (results_directory / f"summary_{mode}.json").write_text(json.dumps(summary, indent=4) + "\n")
    print(json.dumps(summary, indent=2))
    if not isclose(reaction_x, -4000.0, rel_tol=1e-6, abs_tol=1e-3):
        raise RuntimeError(f"{case_directory.name}: support reaction does not balance 4000 N.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("all", *CASES), default="all")
    parser.add_argument("--mode", choices=("average", "center_gp", "center_element"),
                        default="average")
    args = parser.parse_args()
    extractor = load_response_extractor()
    for case in CASES if args.case == "all" else (args.case,):
        print(f"\nRunning {case}, response mode={args.mode}")
        directory = ROOT / case
        model, analysis = create_analysis(directory)
        analysis.Run()
        model_part = model["Structure"]
        values = extractor.extract_von_mises(model_part, args.mode)
        write_results(directory, model_part, values, args.mode)


if __name__ == "__main__":
    main()
