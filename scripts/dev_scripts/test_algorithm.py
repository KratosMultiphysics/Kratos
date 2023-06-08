import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.measurement_likelihood_response_function import MeasurementLikelihoodResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control_system_identification import MaterialPropertiesControlSystemIdentification
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_system_identification import AlgorithmSystemIdentification


os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "measurement_residual_test")))


model = Kratos.Model()
model_part = model.CreateModelPart("Structure")
model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3


# create and fill the optimization problem object
optimization_problem = OptimizationProblem()
optimization_problem.AddProcessType("output_processes")


# creating the execution policy wrapper
execution_policy_wrapper_settings = Kratos.Parameters("""{
            "name"    : "primal_analysis",
            "module": "KratosMultiphysics.OptimizationApplication.execution_policies",
            "type": "stepping_analysis_execution_policy",
            "settings": {
                "model_part_names" : ["Structure"],
                "analysis_module"  : "KratosMultiphysics.StructuralMechanicsApplication",
                "analysis_type"    : "StructuralMechanicsAnalysis",
                "analysis_settings": {
                    "@include_json": "primal_parameters.json"
                }
            },
            "pre_operations"           : [],
            "post_operations"          : [],
            "log_in_file"              : false,
            "log_file_name"            : "structure.log"
        }""")
execution_policy_decorator = ExecutionPolicyDecorator(model, execution_policy_wrapper_settings, optimization_problem)
optimization_problem.AddComponent(execution_policy_decorator)

Kratos.ModelPartIO("model_file", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)


# create control
parameters = Kratos.Parameters("""{
    "model_part_names"      : ["Structure.all_nodes_elements_model_part"],
    "control_variable_name" : "YOUNG_MODULUS",
        "mapping_options":{
            "use_constant_multiplication_factor": true,
            "constant_multiplication_factor" : 10.0
        }
}""")
properties_control = MaterialPropertiesControlSystemIdentification("system_identification_control", model, parameters)
optimization_problem.AddComponent(properties_control)


# create response
response_function_settings = Kratos.Parameters("""{
            "evaluated_model_part_names": ["Structure"],
            "primal_analysis_name"      : "primal_analysis",
            "perturbation_size"         : 1e-8
        }""")
response_function: MeasurementLikelihoodResponseFunction = MeasurementLikelihoodResponseFunction("measurement_residual", model, response_function_settings, optimization_problem)
optimization_problem.AddComponent(response_function)


# Init
execution_policy_decorator.Initialize()
properties_control.Initialize()
response_function.Initialize()

execution_policy_decorator.Execute()


# optimization algorithm object
parameters = Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "AlgorithmSystemIdentification",
            "objective"         : {
                "response_name": "measurement_residual",
                "type"         : "minimization",
                "scaling"      : 1.0
            },
            "controls"          : ["system_identification_control"],
            "echo_level"        : 0,
            "settings"          : {
                "echo_level"      : 0,
                "line_search"     : {
                    "type"          : "const_step",
                    "init_step"     : 0.1,
                    "gradient_scaling": "inf_norm"
                },
                "conv_settings"   : {
                    "type"          : "max_iter",
                    "max_iter"      : 10
                }
            }
        }""")

algorithm = AlgorithmSystemIdentification(model, parameters, optimization_problem)
result = algorithm.SolveOptimizationProblem()
print("Optimization problem solved: {result}")
