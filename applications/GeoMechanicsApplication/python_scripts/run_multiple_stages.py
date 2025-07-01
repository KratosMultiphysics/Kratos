import os
import sys

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers


def run_stages(project_path, n_stages):
    """
    Run all construction stages
    """
    model = Kratos.Model()
    project_parameters_filenames = [f"ProjectParameters_stage{i + 1}.json" for i in range(n_stages)]

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
