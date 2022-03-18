import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.structural_mechanics_analysis_with_co_sim_io import StructuralMechanicsAnalysisWithCoSimIO
from sys import argv

if __name__ == '__main__':
    if len(argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python structural_mechanics_analysis_with_co_sim_io.py <project-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    StructuralMechanicsAnalysisWithCoSimIO(model, parameters).Run()
