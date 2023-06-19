import sys
import os

sys.path.append(os.path.join('..', '..', '..'))

import KratosMultiphysics as Kratos
from KratosMultiphysics.GeoMechanicsApplication import *
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis


def run_stages(project_path,n_stages):
    """
    Run all construction stages

    :param project_path:
    :param n_stages:
    :return:
    """
    parameter_file_names = [os.path.join(project_path, 'ProjectParameters_stage' + str(i + 1) + '.json') for i in
                            range(n_stages)]

    # set stage parameters
    parameters_stages = [None] * n_stages
    os.chdir(project_path)
    for idx, parameter_file_name in enumerate(parameter_file_names):
        with open(parameter_file_name, 'r') as parameter_file:
            parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    stages = [analysis.GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

    [stage.Run() for stage in stages]


if __name__ == "__main__":

    n_stages = int(sys.argv[2])
    project_path = sys.argv[1]

    run_stages(project_path, n_stages)