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
            "youngs_modules"  : [],
            "densities"       : [],
            "beta_settings"   : {
                "initial_value": 25,
                "max_value"    : 25,
                "adaptive"     : false,
                "increase_fac" : 1.5,
                "update_period": 20
            },
            "filter_settings": {
                "type"                      : "explicit",
                "filter_function_type"      : "linear",
                "damping_function_type"     : "sigmoidal",
                "radius"                    : 0.000000000001,
                "max_nodes_in_filter_radius": 1000,
                "linear_solver_settings"    : {},
                "fixed_model_parts"         : []
            }
        }""")
        parameters.RecursivelyValidateAndAssignDefaults(default_settings)

        self.young_modulus = [v for v in parameters["youngs_modules"].GetVector()]
        self.densities =  [v for v in parameters["densities"].GetVector()]
        self.simp_power_fac = parameters["simp_power_fac"].GetInt()
        self.beta = parameters["beta_settings"]["initial_value"].GetDouble()
        self.beta_adaptive = parameters["beta_settings"]["adaptive"].GetBool()
        self.beta_update_period = parameters["beta_settings"]["update_period"].GetInt()
        self.beta_increase_fac = parameters["beta_settings"]["increase_fac"].GetDouble()
        self.beta_max_value = parameters["beta_settings"]["max_value"].GetDouble()

        self.controlled_physical_variables = [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SimpControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

        self.x_values = [i for i, _ in enumerate(self.densities)]

    def Initialize(self) -> None:
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        # init model parts
        self.model_part = self.model_part_operation.GetModelPart()

        # Creating element specific properties for Youngs Modulus and Density
        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        # calculate phi from existing density
        self.phi = self.GetEmptyField()
        KratosOA.PropertiesVariableExpressionIO.Read(self.phi, Kratos.DENSITY)
        self.phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(self.phi, self.x_values, self.densities, self.beta, 1)
        self.__ApplyDensityAndYoungsModulus()

        self.__current_step_value = -1

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

        # now calculate the total sensitivities of density w.r.t. phi
        d_density_d_phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(self.phi, self.x_values, self.densities, self.beta, 1)

        # now calculate the total sensitivities of young modulus w.r.t. phi
        d_young_modulus_d_phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(self.phi, self.x_values, self.young_modulus, self.beta, self.simp_power_fac)

        # now compute response function total sensitivity w.r.t. phi
        d_j_d_phi = d_j_d_density * d_density_d_phi + d_j_d_youngs * d_young_modulus_d_phi
        return d_j_d_phi.Clone()

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        if Kratos.Expression.Utils.NormL2(self.phi - control_field) > 1e-15:
            self.phi = control_field.Clone()
            self.__ApplyDensityAndYoungsModulus()

            step = self.optimization_problem.GetStep()
            # do this only once per step
            if self.__current_step_value != step and self.beta_adaptive:
                self.__current_step_value = step
                if step % self.beta_update_period == 0 and self.beta < self.beta_max_value:
                    self.beta *= self.beta_increase_fac
                    if self.beta > self.beta_max_value:
                        self.beta = self.beta_max_value

    def __ApplyDensityAndYoungsModulus(self) -> None:
        density = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(self.phi, self.x_values, self.densities, self.beta, 1)
        KratosOA.PropertiesVariableExpressionIO.Write(density, Kratos.DENSITY)
        self.un_buffered_data.SetValue("DENSITY", density.Clone(), overwrite=True)

        youngs_modulus = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(self.phi, self.x_values, self.young_modulus, self.beta, self.simp_power_fac)
        KratosOA.PropertiesVariableExpressionIO.Write(youngs_modulus, Kratos.YOUNG_MODULUS)
        self.un_buffered_data.SetValue("YOUNG_MODULUS", youngs_modulus.Clone(), overwrite=True)
