import os
from glob import glob

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.measurement_likelihood_response_function import MeasurementLikelihoodResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control_system_identification import MaterialPropertiesControlSystemIdentification
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_system_identification import AlgorithmSystemIdentification
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis


os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "measurement_residual_test")))

# remove primal analysis files
for file in glob('./Structure_*.vtu'):
    os.remove(file)

with open("optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())
model = Kratos.Model()
analysis = OptimizationAnalysis(model, parameters)
analysis.Run()

# remove primal analysis files
for file in glob('./Structure*.h5'):
    os.remove(file)
