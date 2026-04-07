"""
Mesh convergence analysis script.

Runs simulations for multiple mesh sizes and a reference mesh, then
computes L2 error norms as a post-processing step.

Usage:
    python3.10 mesh_convergence.py
"""
import json
import os

from AnalysisStageWithFlush import CreateAnalysisStageWithFlushInstance
import KratosMultiphysics
import importlib

import scripts.utils as utils


def change_h_value(h_value):
    """Update the mesh size h in config.json."""
    config = utils.get_config()
    config["h"] = h_value
    with open(utils.get_case_dir("config.json"), 'w') as parameter_file:
        json.dump(config, parameter_file, indent=4)


def run_single_simulation(h_value, params_file="ProjectParametersBoussinesq.json"):
    """Run a Kratos simulation for a given mesh size h.

    Updates config.json, generates the mesh if needed, and runs the solver.
    The AnalysisStageWithFlush class handles exporting the solution data
    to HDF5 and MDPA files upon finalization.
    """
    print(f"\n{'='*60}")
    print(f"  Running simulation for h = {h_value}")
    print(f"{'='*60}\n")

    # Update config with the new h value
    change_h_value(h_value)

    with open(params_file, 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

        utils.select_mesh_auto()
        utils.change_peclet(parameters)
        utils.change_output_name_in_params(parameters)

    # Resolve the analysis stage class
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    # Run the simulation
    global_model = KratosMultiphysics.Model()
    simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
    simulation.Run()

    print(f"\n  Simulation for h = {h_value} completed.")


if __name__ == "__main__":
    params_file = "ProjectParametersBoussinesq.json"

    # ── Convergence study parameters ──
    h_values = [0.4, 0.2, 0.1]
    h_ref = 0.05
    # h_values = [0.4, 0.2, 0.1, 0.05]
    # h_ref = 0.025

    results_dir = os.path.join(utils.get_case_dir(), "convergence_results")
    # Remove previous results
    if os.path.exists(results_dir):
        import shutil
        shutil.rmtree(results_dir)

    # ── Phase 1: Run all simulations ──
    all_h = h_values + [h_ref]
    for h_value in all_h:
        run_single_simulation(h_value, params_file)

    # ── Phase 2: Post-processing — compute L2 errors ──
    print(f"\n{'='*60}")
    print("  POST-PROCESSING: Mesh Convergence Analysis")
    print(f"{'='*60}\n")

    from scripts.convergence_postprocess import run_convergence_analysis, plot_convergence

    results = run_convergence_analysis(h_values, h_ref, results_dir)

    # Generate convergence plot
    plot_path = os.path.join(results_dir, "convergence_plot.png")
    plot_convergence(results, output_path=plot_path)

    print(f"\nAll done. Results are in: {results_dir}")
