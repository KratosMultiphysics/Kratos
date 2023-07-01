import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess


with kratos_unittest.WorkFolderScope(".", __file__):
    with open("optimization_parameters.json", "r") as file_input:
        parameters = Kratos.Parameters(file_input.read())
    model = Kratos.Model()
    analysis = OptimizationAnalysis(model, parameters)
    analysis.Run()
