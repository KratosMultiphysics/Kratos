import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"ShellThicknessControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"ShellThicknessControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return ShellThicknessControl(parameters["name"].GetString(), model, parameters["settings"])

def GetImplicitFilterParameters(model_part: Kratos.ModelPart, parameters: Kratos.Parameters):
    model_part_name = model_part.FullName()
    filter_radius = parameters["filter_settings"]["radius"].GetDouble()
    linear_solver_settings = parameters["filter_settings"]["linear_solver_settings"]
    shell_scalar_filter_parameters = Kratos.Parameters("""
     {
        "solver_settings" : {
             "domain_size"          : 3,
             "echo_level"           : 0,
             "filter_type"          : "general_scalar",
             "model_import_settings": {
                 "input_type"     : "use_input_model_part"
             }
         },
        "problem_data": {
             "echo_level"    : 0,
             "start_time"    : 0.0,
             "end_time"      : 1.0,
             "parallel_type" : "OpenMP"
         },
        "processes" : {
                "boundary_conditions_process_list" : []
            }
     }""")

    shell_scalar_filter_parameters["solver_settings"].AddDouble("filter_radius",filter_radius)
    shell_scalar_filter_parameters["solver_settings"].AddString("model_part_name",model_part_name)
    shell_scalar_filter_parameters["solver_settings"].AddValue("linear_solver_settings", linear_solver_settings)

    fixed_model_parts = parameters["fixed_model_parts"].GetStringArray()

    for fixed_model_part in fixed_model_parts:
        auto_process_settings = Kratos.Parameters(
            """
            {
                "python_module" : "fix_scalar_variable_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "FixScalarVariableProcess",
                "Parameters"    : {
                    "model_part_name" : \""""+fixed_model_part+"""\",
                    "variable_name"   : "HELMHOLTZ_SCALAR",
                    "constrained"     : true
                }
            }
            """)
        shell_scalar_filter_parameters["processes"]["boundary_conditions_process_list"].Append(auto_process_settings)

    return shell_scalar_filter_parameters

class ShellThicknessControl(Control):
    """Shell thickness control

    This is filtering-based discrete thickness control which parametrizes and controls
    shell elements thickness.

    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        self.parameters = parameters
        self.model = model

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [],
            "filter_settings": {
                "type" : "implicit",
                "radius": 0.000000000001,
                "linear_solver_settings" : {}
            },
            "beta_settings": {
                "initial_value" : 25,
                "max_value" : 25,
                "adaptive" : false,
                "increase_fac" : 1.5,
                "update_period" : 20
            },
            "SIMP_power_fac": 3,
            "initial_thickness":0.000001,
            "physical_thicknesses": [],
            "fixed_model_parts": [],
            "fixed_model_parts_thicknesses": [],
            "utilities": []
        }""")

        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.parameters["beta_settings"].ValidateAndAssignDefaults(default_settings["beta_settings"])
        self.parameters["filter_settings"].ValidateAndAssignDefaults(default_settings["filter_settings"])

        self.filter_type = self.parameters["filter_settings"]["type"].GetString()
        self.supported_filter_types = ["implicit"]

        if not self.filter_type in self.supported_filter_types:
            raise RuntimeError(f"The specified filter type \"{self.filter_type}\" is not supported. Followings are the variables:\n\t" + "\n\t".join([k for k in self.supported_filter_types]))

        controlled_model_parts = [model[model_part_name] for model_part_name in parameters["controlled_model_part_names"].GetStringArray()]
        if len(controlled_model_parts) == 0:
            raise RuntimeError(f"No model parts were provided for ShellThicknessControl. [ control name = \"{self.GetName()}\"]")

        root_model_part = controlled_model_parts[0].GetRootModelPart()
        self.model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.UNION, root_model_part, controlled_model_parts, False)

        self.fixed_model_parts = self.parameters["fixed_model_parts"].GetStringArray()
        self.fixed_model_parts_thicknesses = self.parameters["fixed_model_parts_thicknesses"].GetVector()
        if len(self.fixed_model_parts) != len(self.fixed_model_parts_thicknesses):
            raise RuntimeError(f"fixed_model_parts and fixed_model_parts_thicknesses should have the same size !")

        helmholtz_settings = GetImplicitFilterParameters(self.model_part, self.parameters)

        self.filter = HelmholtzAnalysis(self.model, helmholtz_settings)

    def Initialize(self) -> None:
        ModelPartUtilities.ExecuteOperationOnModelPart(self.model_part)

        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        self.filter.Initialize()

        # create the control field
        self.control_field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(self.control_field, self.parameters["initial_thickness"].GetDouble())

        self._SetFixedModelPartValues()
        thickness_physical_field = self.filter.FilterField(self.control_field)

        # now update physical field
        KratosOA.PropertiesVariableExpressionIO.Write(thickness_physical_field, Kratos.THICKNESS)

    def Check(self):
        return self.filter.Check()

    def Finalize(self):
        return self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.THICKNESS]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        with TimeLogger(self.__class__.__name__, f"Mapping Gardient {self.GetName()}...", f"Finished updating of {self.GetName()}."):
            keys = physical_gradient_variable_container_expression_map.keys()
            if len(keys) != 1:
                raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
            if Kratos.THICKNESS not in keys:
                raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {Kratos.THICKNESS.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

            physical_gradient = physical_gradient_variable_container_expression_map[Kratos.THICKNESS]
            if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
                raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

            physical_gradient_variable_container_expression_map[Kratos.THICKNESS].Clone()
            self._SetFixedModelPartValues(True)
            filtered_gradient = self.filter.FilterIntegratedField(physical_gradient_variable_container_expression_map[Kratos.THICKNESS])

            return filtered_gradient

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")

        if KratosOA.ExpressionUtils.NormL2(self.control_field - new_control_field) > 1e-9:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}."):
                self._SetFixedModelPartValues()
                new_physical_field = self.filter.FilterField(new_control_field)

                # now update physical field
                KratosOA.PropertiesVariableExpressionIO.Write(new_physical_field, Kratos.THICKNESS)
                return True

        return False

    def _SetFixedModelPartValues(self, is_backward = False) -> None:
        for model_part_name, model_part_value in zip(self.fixed_model_parts, self.fixed_model_parts_thicknesses):
            if is_backward:
                model_part_value = 0.0
            field =  Kratos.Expression.NodalExpression(self.model[model_part_name])
            Kratos.Expression.LiteralExpressionIO.SetData(field, model_part_value)
            Kratos.Expression.VariableExpressionIO.Write(field, KratosOA.HELMHOLTZ_SCALAR, True)

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {Kratos.THICKNESS.Name()}"
