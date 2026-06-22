import sys
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers


def _validated_stage_file_paths(input_path, n_stages, filename_pattern, base_path=None):
    """Build stage file paths and ensure they stay within the configured input directory."""
    base_path = Path(base_path).expanduser().resolve() if base_path is not None else None
    input_dir = Path(input_path).expanduser()
    if not input_dir.is_absolute() and base_path is not None:
        input_dir = base_path / input_dir
    input_dir = input_dir.resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input path does not exist or is not a directory: {input_dir}")

    result = []
    for i in range(n_stages):
        candidate = (input_dir / filename_pattern.format(i + 1)).resolve()
        try:
            candidate.relative_to(input_dir)
        except ValueError as exc:
            raise ValueError(
                f"Stage file path escapes input path boundary: {candidate}"
            ) from exc
        result.append(candidate)

    return result


def run_stages(
    project_path,
    n_stages,
    filename_pattern="ProjectParameters_stage{}.json",
    input_path=None,
):
    """
    Run all construction stages
    """
    model = Kratos.Model()

    if input_path is None:
        input_path = project_path

    project_parameters_filenames = _validated_stage_file_paths(
        input_path, n_stages, filename_pattern, base_path=project_path
    )

    result = []
    with context_managers.set_cwd_to(project_path):
        for filename in project_parameters_filenames:
            with open(filename, "r") as f:
                stage_parameters = Kratos.Parameters(f.read())
            result.append(analysis.GeoMechanicsAnalysis(model, stage_parameters))
            result[-1].Run()

    return result


if __name__ == "__main__":
    n_stages = int(sys.argv[2])
    project_path = sys.argv[1]

    run_stages(project_path, n_stages)
