import KratosMultiphysics

from KratosMultiphysics.LaserDrillingApplication.laserdrilling_analysis import LaserDrillingAnalysis

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg = "Too many input arguments!\n"
        err_msg += "Use this script in the following way:\n"
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 laserdrilling_analysis.py"\n'
        err_msg += "- With custom ProjectParameters:\n"
        err_msg += '    "python3 laserdrilling_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2:  # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else:  # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name, "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = LaserDrillingAnalysis(model, parameters)
    simulation.Run()
