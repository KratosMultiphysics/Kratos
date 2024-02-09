import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"SimpControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SimpControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SimpControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class SimpControl(Control):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [""],
            "simp_power_fac"  : 3,
            "filter_settings": {
                "type"                      : "explicit",
                "filter_function_type"      : "linear",
                "damping_function_type"     : "sigmoidal",
                "radius"                    : 0.000000000001,
                "max_nodes_in_filter_radius": 1000,
                "fixed_model_parts"         : []
            }
        }""")
        parameters.RecursivelyValidateAndAssignDefaults(default_settings)
        self.simp_power_fac = parameters["simp_power_fac"].GetDouble()

        self.controlled_physical_variables = [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

        # filtering settings
        filter_settings = parameters["filter_settings"]
        filter_type = filter_settings["type"].GetString()
        if filter_type != "explicit":
            raise RuntimeError(f"Currently simp control only support \"explicit\" vertex morphing.")

        self.filter_function_type = filter_settings["filter_function_type"].GetString()
        self.damping_function_type = filter_settings["damping_function_type"].GetString()
        self.filter_radius = filter_settings["radius"].GetDouble()
        self.max_nodes_in_filter_radius = filter_settings["max_nodes_in_filter_radius"].GetInt()
        self.fixed_model_parts = filter_settings["fixed_model_parts"].GetStringArray()

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SimpControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

    def Initialize(self) -> None:
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        # init model parts
        self.model_part = self.model_part_operation.GetModelPart()

        # Creating element specific properties for Youngs Modulus and Density
        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        # calculate phi from existing density
        self.phi = Kratos.Expression.Utils.Collapse(self.GetEmptyField() + 1.0)
        self.density = 0.0
        self.young_modulus = 0.0
        for element in self.model_part.Elements:
            self.density = element.Properties[Kratos.DENSITY]
            self.young_modulus = element.Properties[Kratos.YOUNG_MODULUS]
            break

        self.__ApplyDensityAndYoungsModulus()

        # create vertex morphing filter
        if self.fixed_model_parts:
            fixed_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", self.fixed_model_parts, False)
            fixed_model_part = fixed_model_part_operation.GetModelPart()
            self.filter = KratosOA.ElementExplicitFilter(self.model_part, fixed_model_part, self.filter_function_type, self.damping_function_type, self.max_nodes_in_filter_radius)
            filter_radius = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, self.filter_radius)
            self.filter.SetFilterRadius(filter_radius)
        else:
            self.filter = KratosOA.ElementExplicitFilter(self.model_part, self.filter_function_type, self.max_nodes_in_filter_radius)
            filter_radius = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, self.filter_radius)
            self.filter.SetFilterRadius(filter_radius)

        self.un_buffered_data.SetValue("filtered_phi", self.phi, overwrite=True)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.controlled_physical_variables

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.phi

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 2:
            raise RuntimeError(f"Not provided required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if  Kratos.DENSITY not in keys or Kratos.YOUNG_MODULUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. DENSITY or YOUNG_MODULUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the density partial sensitivity of the response function
        d_j_d_density = physical_gradient_variable_container_expression_map[Kratos.DENSITY]

        # second calculate the young modulus partial sensitivity of the response function
        d_j_d_youngs = physical_gradient_variable_container_expression_map[Kratos.YOUNG_MODULUS]

        # now calculate the total sensitivities of young modulus w.r.t. phi
        d_young_modulus_d_phi = Kratos.Expression.Utils.Pow(self.phi, self.simp_power_fac - 1) * self.young_modulus * self.simp_power_fac

        # now compute response function total sensitivity w.r.t. phi
        d_j_d_phi = d_j_d_density * self.density + d_j_d_youngs * d_young_modulus_d_phi

        return self.filter.FilterIntegratedField(d_j_d_phi)

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        if Kratos.Expression.Utils.NormL2(self.phi - control_field) > 1e-15:
            self.phi = KratosOA.ExpressionUtils.Clamp(self.filter.FilterField(control_field), 0, 1)
            self.un_buffered_data.SetValue("filtered_phi", self.phi.Clone(), overwrite=True)
            self.__ApplyDensityAndYoungsModulus()

    def __ApplyDensityAndYoungsModulus(self) -> None:
        density = self.phi * self.density
        KratosOA.PropertiesVariableExpressionIO.Write(density, Kratos.DENSITY)
        self.un_buffered_data.SetValue("DENSITY", density.Clone(), overwrite=True)

        youngs_modulus = Kratos.Expression.Utils.Pow(self.phi, self.simp_power_fac) * self.young_modulus
        KratosOA.PropertiesVariableExpressionIO.Write(youngs_modulus, Kratos.YOUNG_MODULUS)
        self.un_buffered_data.SetValue("YOUNG_MODULUS", youngs_modulus.Clone(), overwrite=True)
