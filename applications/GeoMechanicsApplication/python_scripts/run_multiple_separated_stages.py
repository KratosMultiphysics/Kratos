import sys
import os
import KratosMultiphysics

from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

if __name__ == "__main__":

    model = KratosMultiphysics.Model()
    project_path = sys.argv[1]
    n_stages = int(sys.argv[2])

    directory_names = [os.path.join(project_path, 'Stage_' +  str(i + 1)) for i in range(n_stages)]

    # set stage parameters
    parameters_stages = [None] * n_stages
    for idx, directory_name in enumerate(directory_names):
        parameter_file_name = directory_name + '/ProjectParameters.json'
        with open(parameter_file_name, 'r') as parameter_file:
            parameters_stages[idx] = KratosMultiphysics.Parameters(parameter_file.read())

    stages = [GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

    for idx, stage in enumerate(stages):
        os.chdir(directory_names[idx])
        stage.Run()

