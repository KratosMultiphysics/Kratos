import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.measurement_likelihood_response_function import MeasurementLikelihoodResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.material.material_properties_control_system_identification import MaterialPropertiesControlSystemIdentification

model = Kratos.Model()
model_part = model.CreateModelPart("Structure")
model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
with kratos_unittest.WorkFolderScope("measurement_residual_test", __file__):
    Kratos.ModelPartIO("model_file", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)
    material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "StructuralMaterials.json"}} """)
    Kratos.ReadMaterialsUtility(material_settings, model)
parameters = Kratos.Parameters("""{
    "model_part_names"      : ["Structure.all_nodes_elements_model_part"],
    "control_variable_name" : "YOUNG_MODULUS"
}""")
properties_control = MaterialPropertiesControlSystemIdentification("system_identification_control", model, parameters)

properties_control.Initialize()

for element in model_part.Elements:
    print(f"Previous: {element.Properties[Kratos.YOUNG_MODULUS]}")
    element.Properties[Kratos.YOUNG_MODULUS] += 12

control_model_part = model_part.GetSubModelPart("Union_Structure#all_nodes_elements_model_part_EN")

temp = KratosOA.ContainerExpression.ElementPropertiesExpression(control_model_part)
temp.Read(Kratos.YOUNG_MODULUS)
properties_control.Update(temp)

for element in model_part.Elements:
    print(f"Updated: {element.Properties[Kratos.YOUNG_MODULUS]}")

field_to_be_mapped = {Kratos.YOUNG_MODULUS: KratosOA.ContainerExpression.ElementPropertiesExpression(control_model_part)}
mapped_container = properties_control.MapGradient(field_to_be_mapped)

for element in model_part.Elements:
    print(f"Mapped: {element.Properties[Kratos.YOUNG_MODULUS]}")
