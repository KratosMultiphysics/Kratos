import os
import KratosMultiphysics as Kratos

from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

if __name__ == "__main__":

    currentWorking = os.getcwd()

    # construct parameterfile names of stages to run
    project_path = r"C:\Users\graaf\OneDrive - Stichting Deltares\Documents\settlement_app\scenario_1\input_files"
    n_stages = 4
    parameter_file_names = [os.path.join(project_path, 'ProjectParameters_stage' + str(i + 1) + '.json') for i in
                            range(n_stages)]

    # change to project directory
    os.chdir(project_path)

    # setup stages from parameterfiles
    parameters_stages = [None] * n_stages
    for idx, parameter_file_name in enumerate(parameter_file_names):
        with open(parameter_file_name, 'r') as parameter_file:
            parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    stages = [GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

    # execute the stages
    for stage in stages:
        stage.Run()

    # back to working directory
    os.chdir(currentWorking)
