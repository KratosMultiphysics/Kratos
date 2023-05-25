import json
import dataclasses
from typing import List

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.measurement_likelihood_response_function import MeasurementLikelihoodResponseFunction

model = Kratos.Model()
optimization_problem = OptimizationProblem()
model_part = model.CreateModelPart("Structure")
model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

# create the primal analysis execution policy wrapper
with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
        # creating the execution policy wrapper
        execution_policy_wrapper_settings = Kratos.Parameters("""{
            "name"    : "primal",
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

        Kratos.ModelPartIO("Structure", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)

        # creating the response function wrapper
        response_function_settings = Kratos.Parameters("""{
            "evaluated_model_part_names": ["Structure"],
            "primal_analysis_name"      : "primal",
            "perturbation_size"         : 1e-8
        }""")
        response_function: MeasurementLikelihoodResponseFunction = MeasurementLikelihoodResponseFunction("strain_energy", model, response_function_settings, optimization_problem)
        optimization_problem.AddComponent(response_function)

        execution_policy_decorator.Initialize()
        response_function.Initialize()

        # now replace the properties
        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model["Structure.structure"], model_part.Elements)

        execution_policy_decorator.Execute()
        ref_value = response_function.CalculateValue()
        print(f"Objective value {ref_value}")

        sensitivity_expression = KratosOA.ContainerExpression.CollectiveExpressions([KratosOA.ContainerExpression.ElementPropertiesExpression(model_part)])
        response_function.CalculateGradient({Kratos.YOUNG_MODULUS: sensitivity_expression})
        for element in model_part.Elements:
                print(f"Sensitivity: {element.Properties[KratosOA.YOUNG_MODULUS_SENSITIVITY]}")

        # file = open("./measurements.json")
        # data = json.load(file)
        # print(data)
@dataclasses.dataclass
class LoadDataContainer:
        type_of_load:str = "PointLoad",
        direction_normal:List[float] = [0,-1,0],
        strength_in_N:float = 1,
        position_of_mesh_vertex:List[float] = [0,0,0],
@dataclasses.dataclass
class SensorDataContainer:
        type_of_sensor:str = "Displacement",
        measurement_direction_normal:List[float] = [0,1,0],
        position_of_mesh_vertex:List[float] = [0,0,0],
        measured_value:float = 0
@dataclasses.dataclass
class MeasurementDataContainer:
        description_of_load:LoadDataContainer = None,
        description_of_sensors:List[SensorDataContainer] = None,

load = LoadDataContainer(strength_in_N=2000,
                         position_of_mesh_vertex=[0,0,0])

sensor = SensorDataContainer(measured_value=0.5,
                             position_of_mesh_vertex=[0,0,0])

data = MeasurementDataContainer(description_of_load=load,
                                description_of_sensors=[sensor])

with open(".data.json", "w") as outfile:
    json.dump(dataclasses.asdict(data), outfile)

file = open(".data.json")
data2 = json.load(file)
print(data2)
