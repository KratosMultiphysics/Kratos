import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

if __name__ == "__main__":
    model = KratosMultiphysics.Model()
    with open("auxiliary_files/system_identification_p_norm/optimization_parameters.json", "r") as file_input:
        params = KratosMultiphysics.Parameters(file_input.read())

    analysis = OptimizationAnalysis(model, params)
    analysis.Run()