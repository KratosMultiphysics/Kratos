import sys
import os

sys.path.append(os.path.join("..", "..", ".."))

import KratosMultiphysics as Kratos
from KratosMultiphysics.GeoMechanicsApplication import *
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis


def run_stages(project_path, n_stages):
    """
    Run all construction stages
    """
    os.chdir(project_path)

    all_stage_parameters = []
    for i in range(n_stages):
        with open(f"ProjectParameters_stage{i + 1}.json", "r") as f:
            all_stage_parameters.append(Kratos.Parameters(f.read()))

    model = Kratos.Model()

    for stage_parameters in all_stage_parameters:
        stage = analysis.GeoMechanicsAnalysis(model, stage_parameters)
        stage.Run()


if __name__ == "__main__":
    n_stages = int(sys.argv[2])
    project_path = sys.argv[1]

    run_stages(project_path, n_stages)
