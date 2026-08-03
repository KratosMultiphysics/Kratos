import sys
import os
from pathlib import Path

import KratosMultiphysics
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

def _validated_stage_directory_paths(project_path, n_stages):
    """
    Builds absolute paths for each stage directory and ensures they stay 
    within the configured project directory.
    """
    # Define the immutable sandbox based on where the user pointed the CLI
    safe_base = Path(project_path).expanduser().resolve()
    if not safe_base.is_dir():
        raise FileNotFoundError(f"Project path does not exist or is not a directory: {safe_base}")

    if(n_stages <=0):
        raise ValueError(f"Number of stages must be positive: {n_stages}")
    
    result = []
    for i in range(n_stages):
        # Construct the candidate directory name
        dir_name = f"Stage_{i + 1}"
        candidate = (safe_base / dir_name).resolve()
        # Security Check 1: Ensure the resolved stage directory stays within the project root (e.g., prevent symlink escape)
        try:
            candidate.relative_to(safe_base)
        except ValueError as exc:
            raise ValueError(
                f"Stage directory '{dir_name}' attempts to escape the project root."
            ) from exc
        # Security Check 2: Ensure the physical directory actually exists
        if not candidate.is_dir():
            raise FileNotFoundError(f"Required stage directory not found: {candidate}")
            
        result.append(candidate)

    return result

if __name__ == "__main__":

    model = KratosMultiphysics.Model()
    project_path = sys.argv[1]
    n_stages = int(sys.argv[2])

    stage_directories = _validated_stage_directory_paths(project_path, n_stages)

    parameters_stages = [
        KratosMultiphysics.Parameters(
           (stage_directory / "ProjectParameters.json").read_text(encoding='utf-8')
        )
        for stage_directory in stage_directories
    ]

    stages = [GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

    for idx, stage in enumerate(stages):
        os.chdir(stage_directories[idx])
        stage.Run()

