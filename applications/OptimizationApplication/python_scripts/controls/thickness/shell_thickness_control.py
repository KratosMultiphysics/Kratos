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
    if not parameters.Has("controlled_model_part_names"):
        raise RuntimeError(f"ShellThicknessControl instantiation requires a \"controlled_model_part_names\" in parameters [ parameters = {parameters}].")

    return ShellThicknessControl(parameters["name"].GetString(), model, parameters)

def GetImplicitFilterParameters(parameters):
    filter_radius = parameters["settings"]["filter_settings"]["radius"].GetDouble()
    model_part_name = parameters["controlled_model_part_names"].GetStringArray()[0]
    linear_solver_settings = parameters["settings"]["filter_settings"]["linear_solver_settings"]
    shell_scalar_filter_parameters = Kratos.Parameters("""
     {
        "solver_settings" : {
             "domain_size"     : 3,
             "echo_level"      : 0,
             "filter_type"     : "general_scalar",
             "model_import_settings"              : {
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

    fixed_model_parts = parameters["settings"]["fixed_model_parts"].GetStringArray()

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
            "filter_settings": {},
            "beta_settings": {},
            "SIMP_power_fac": 3,
            "initial_thickness":0.000001,
            "physical_thicknesses": [],
            "fixed_model_parts": [],
            "fixed_model_parts_thicknesses": [],
            "utilities": []
        }""")

        self.parameters["settings"].ValidateAndAssignDefaults(default_settings)

        beta_default_setting = Kratos.Parameters("""{
            "initial_value" : 25,
            "max_value" : 25,
            "adaptive" : false,
            "increase_fac" : 1.5,
            "update_period" : 20
        }""")

        self.parameters["settings"]["beta_settings"].ValidateAndAssignDefaults(beta_default_setting)

        filter_default_setting = Kratos.Parameters("""{
            "type" : "implicit",
            "radius": 0.000000000001,
            "linear_solver_settings" : {}
        }""")
        self.parameters["settings"]["filter_settings"].ValidateAndAssignDefaults(filter_default_setting)
        self.filter_type = self.parameters["settings"]["filter_settings"]["type"].GetString()
        self.supported_filter_types = ["implicit"]

        if not self.filter_type in self.supported_filter_types:
            raise RuntimeError(f"The specified filter type \"{self.filter_type}\" is not supported. Followings are the variables:\n\t" + "\n\t".join([k for k in self.supported_filter_types]))

        linear_solver_default_setting = Kratos.Parameters("""{
            "solver_type" : "amgcl",
            "smoother_type":"ilu0",
            "krylov_type": "gmres",
            "coarsening_type": "aggregation",
            "max_iteration": 200,
            "provide_coordinates": false,
            "gmres_krylov_space_dimension": 100,
            "verbosity" : 0,
            "tolerance": 1e-7,
            "scaling": false,
            "block_size": 1,
            "use_block_matrices_if_possible" : true,
            "coarse_enough" : 5000
        }""")

        self.parameters["settings"]["filter_settings"]["linear_solver_settings"].ValidateAndAssignDefaults(linear_solver_default_setting)

        self.model_part_names = parameters["controlled_model_part_names"].GetStringArray()
        # TODO : we need to modify UnionModelParts such that it does not give error
        # remove this part after modifying UnionModelParts
        if len(self.model_part_names)>1:
            raise RuntimeError(f"A single model part is allowed")

        self.model_part_name = self.model_part_names[0]
        self.model_part = self.model.GetModelPart(self.model_part_name)

        self.fixed_model_parts = self.parameters["settings"]["fixed_model_parts"].GetStringArray()
        self.fixed_model_parts_thicknesses = self.parameters["settings"]["fixed_model_parts_thicknesses"].GetVector()
        if len(self.fixed_model_parts) != len(self.fixed_model_parts_thicknesses):
            raise RuntimeError(f"fixed_model_parts and fixed_model_parts_thicknesses should have the same size !")

        helmholtz_settings = GetImplicitFilterParameters(self.parameters)

        self.filter = HelmholtzAnalysis(self.model, helmholtz_settings)

        self.controlled_physical_variable = Kratos.KratosGlobals.GetVariable("THICKNESS")

    def Initialize(self) -> None:

        # TODO. after solving issue with UnionModelParts such that it creates the union model part eventhough it exists
        # we uncomment the following lines and work with the unioned nodel part
        # model_parts_list = [self.model[model_part_name] for model_part_name in self.model_part_names]
        # root_model_part = model_parts_list[0].GetRootModelPart()
        # is_new_model_part, self.model_part = ModelPartUtilities.UnionModelParts(root_model_part, model_parts_list, False)

        # if is_new_model_part:
            # now create entity specific properties for the merged model part which is used for the control.
        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)

        self.filter.Initialize()

        # create the control field
        self.control_field = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)
        self.control_field.SetData(self.parameters["settings"]["initial_thickness"].GetDouble())
        self.SetFixedModelPartValues()
        thickness_physical_field = self.filter.FilterField(self.control_field)
        # now update physical field
        thickness_physical_field.Evaluate(self.controlled_physical_variable)
        return True

    def Check(self):
        return self.filter.Check()

    def Finalize(self):
        return self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [self.controlled_physical_variable]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)
        field.SetData(0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field

    def MapGradient(self, physical_gradient_variable_container_expression_map: dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]) -> ContainerExpressionTypes:
        with TimeLogger(self.__class__.__name__, f"Mapping Gardient {self.GetName()}...", f"Finished updating of {self.GetName()}."):
            keys = physical_gradient_variable_container_expression_map.keys()
            if len(keys) != 1:
                raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
            if self.controlled_physical_variable not in keys:
                raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {self.controlled_physical_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

            physical_gradient = physical_gradient_variable_container_expression_map[self.controlled_physical_variable]
            if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
                raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

            physical_gradient_variable_container_expression_map[self.controlled_physical_variable].Clone()
            self.SetFixedModelPartValues(True)
            filtered_gradient = self.filter.FilterIntegratedField(physical_gradient_variable_container_expression_map[self.controlled_physical_variable])

            return filtered_gradient

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")
        if KratosOA.ContainerExpressionUtils.NormL2(self.control_field - new_control_field) > 1e-9:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}."):
                self.SetFixedModelPartValues()
                new_physical_field = self.filter.FilterField(new_control_field)
                # now update physical field
                new_physical_field.Evaluate(self.controlled_physical_variable)
                return True

        return False

    def SetFixedModelPartValues(self, is_bachward = False):

        for i in range(len(self.fixed_model_parts)):
            model_part_name = self.fixed_model_parts[i]
            model_part_value = self.fixed_model_parts_thicknesses[i]
            if is_bachward:
                model_part_value = 0.0
            field =  Kratos.ContainerExpression.HistoricalExpression(self.model.GetModelPart(model_part_name))
            field.SetData(model_part_value)
            field.Evaluate(KratosOA.HELMHOLTZ_SCALAR)

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"
