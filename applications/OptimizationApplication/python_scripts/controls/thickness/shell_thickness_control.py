from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"ShellThicknessControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"ShellThicknessControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return ShellThicknessControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

def GetImplicitFilterParameters(model_part_name: str, parameters: Kratos.Parameters) -> Kratos.Parameters:
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

    fixed_model_parts = parameters["fixed_model_parts_and_thicknesses"].keys()

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
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [],
            "filter_settings": {
                "type" : "implicit",
                "radius": 0.000000000001,
                "linear_solver_settings" : {}
            },
            "beta_value": 25.0,
            "penalty_power": 1,
            "output_all_fields": false,
            "initial_physical_thickness":0.000001,
            "physical_thicknesses": [],
            "fixed_model_parts_and_thicknesses": {},
            "utilities": []
        }""")

        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.parameters["filter_settings"].ValidateAndAssignDefaults(default_settings["filter_settings"])

        self.penalty_power = self.parameters["penalty_power"].GetInt()
        self.beta = self.parameters["beta_value"].GetDouble()
        self.output_all_fields = self.parameters["output_all_fields"].GetBool()
        self.physical_thicknesses = self.parameters["physical_thicknesses"].GetVector()
        num_phys_thick = len(self.physical_thicknesses)
        if num_phys_thick == 0:
            raise RuntimeError("The physical_thicknesses can not be empty, at least min and max should be provided.")
        self.filtered_thicknesses = [i for i, _ in enumerate(self.physical_thicknesses)]

        self.filter_type = self.parameters["filter_settings"]["type"].GetString()
        self.supported_filter_types = ["implicit"]

        if not self.filter_type in self.supported_filter_types:
            raise RuntimeError(f"The specified filter type \"{self.filter_type}\" is not supported. Followings are the variables:\n\t" + "\n\t".join([k for k in self.supported_filter_types]))

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts were provided for ShellThicknessControl. [ control name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'Optional[Kratos.ModelPart]' = None

        self.fixed_model_parts_and_thicknesses: 'dict[str, float]' = {}
        for model_part_name, thickness_params in self.parameters["fixed_model_parts_and_thicknesses"].items():
            self.fixed_model_parts_and_thicknesses[model_part_name] = thickness_params.GetDouble()

        if not all(entry.GetDouble() in self.physical_thicknesses for entry in self.fixed_model_parts_and_thicknesses.values()):
            raise RuntimeError("fixed_model_parts_thicknesses should exist in physical_thicknesses !")

        helmholtz_settings = GetImplicitFilterParameters(self.model_part_operation.GetModelPartFullName(), self.parameters)

        # create implicit filter using Helmholtz analysis
        self.filter = HelmholtzAnalysis(self.model, helmholtz_settings)

        self.is_initialized = False

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        # initiliaze the filter
        self.filter.Initialize()
        # compute control thickness field from the initial physical control field
        physical_thickness_field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(physical_thickness_field, self.parameters["initial_physical_thickness"].GetDouble())
        # project backward the uniform physical control field and assign it to the control field
        self.control_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(
                                                physical_thickness_field,
                                                self.filtered_thicknesses,
                                                self.physical_thicknesses,
                                                self.beta,
                                                self.penalty_power)
        # now update the physical thickness field
        self._UpdateAndOutputFields()
        self.is_initialized = True

    def Check(self) -> None:
        return self.filter.Check()

    def Finalize(self) -> None:
        return self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.THICKNESS]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field

    def GetPhysicalField(self) -> ContainerExpressionTypes:
        physical_thickness_field = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(physical_thickness_field, Kratos.THICKNESS)
        return physical_thickness_field

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        with TimeLogger("ShellThicknessControl::MapGradient", None, "Finished",False):
            keys = physical_gradient_variable_container_expression_map.keys()
            if len(keys) != 1:
                raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
            if Kratos.THICKNESS not in keys:
                raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {Kratos.THICKNESS.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

            physical_gradient = physical_gradient_variable_container_expression_map[Kratos.THICKNESS]
            if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
                raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

            # multiply the physical sensitivity field with projection derivatives
            projected_gradient = physical_gradient * self.projection_derivative_field
            # now filter the field
            filtered_gradient = self.filter.FilterIntegratedField(projected_gradient)

            return filtered_gradient

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")

        if Kratos.Expression.Utils.NormL2(self.control_field - new_control_field) > 1e-15:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}.",False):
                # update the control thickness field
                self.control_field = new_control_field
                # now update the physical field
                self._UpdateAndOutputFields()
                return True
        return False

    def _UpdateAndOutputFields(self) -> None:
        # apply boundary conditions for the forward filtering
        self._SetFixedModelPartValues()
        # filter the control field
        filtered_thickness_field = self.filter.FilterField(self.control_field)
        # project forward the filtered thickness field
        projected_filtered_thickness_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(
                                                filtered_thickness_field,
                                                self.filtered_thicknesses,
                                                self.physical_thicknesses,
                                                self.beta,
                                                self.penalty_power)
        # now update physical field
        KratosOA.PropertiesVariableExpressionIO.Write(projected_filtered_thickness_field, Kratos.THICKNESS)
        # compute and strore projection derivatives for consistent filtering of the sensitivities
        self.projection_derivative_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(
                                                filtered_thickness_field,
                                                self.filtered_thicknesses,
                                                self.physical_thicknesses,
                                                self.beta,
                                                self.penalty_power)
        # apply boundary conditions for the filtering of the sensitivites
        self._SetFixedModelPartValues(True)
        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("physical_thickness", projected_filtered_thickness_field.Clone(),overwrite=self.is_initialized)
        if self.output_all_fields:
            un_buffered_data.SetValue("filtered_thickness", filtered_thickness_field.Clone(),overwrite=self.is_initialized)
            un_buffered_data.SetValue("control_thickness", self.control_field.Clone(),overwrite=self.is_initialized)
            un_buffered_data.SetValue("projection_derivative", self.projection_derivative_field.Clone(),overwrite=self.is_initialized)

    def _SetFixedModelPartValues(self, is_backward = False) -> None:
        for model_parts_name, model_parts_thickness in self.fixed_model_parts_and_thicknesses.items():
            if is_backward:
                model_parts_thickness = 0.0
            field =  Kratos.Expression.NodalExpression(self.model[model_parts_name])
            Kratos.Expression.LiteralExpressionIO.SetData(field, model_parts_thickness)
            Kratos.Expression.VariableExpressionIO.Write(field, KratosOA.HELMHOLTZ_SCALAR, True)

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {Kratos.THICKNESS.Name()}"
