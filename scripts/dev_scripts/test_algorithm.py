import os
import shutil
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


def pre_operations():
    try:
        shutil.rmtree("./vtu_results")
    except:
        pass
    try:
        shutil.rmtree("./other_results")
    except:
        pass


def post_operations():
    os.makedirs("./vtu_results")
    os.makedirs("./other_results")

    # remove primal analysis files
    for file in glob('./Structure*.h5'):
        os.remove(file)

    # move results to specific folders
    for file in glob('./Structure_*.vtu'):
        os.rename(str(file), "./vtu_results"+str(file)[1:])

    for file in glob('./*.csv')+glob('./*.time')+glob('./*.html'):
        os.rename(str(file), "./other_results"+str(file)[1:])


working_folder = "measurement_residual_test"
os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

pre_operations()

with open("optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())
model = Kratos.Model()
analysis = OptimizationAnalysis(model, parameters)
analysis.Run()

post_operations()
