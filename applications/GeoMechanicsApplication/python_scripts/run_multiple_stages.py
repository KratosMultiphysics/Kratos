import os
import sys

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers


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

    project_parameters_filenames = [
        os.path.join(input_path, filename_pattern.format(i + 1))
        for i in range(n_stages)
    ]

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
