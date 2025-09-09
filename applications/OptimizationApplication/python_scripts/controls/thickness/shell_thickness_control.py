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
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.opt_projection import CreateProjection

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"ShellThicknessControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"ShellThicknessControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return ShellThicknessControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

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
            "controlled_model_part_names"       : [],
            "filter_settings"                   : {},
            "output_all_fields"                 : false,
            "physical_thicknesses"              : [],
            "consider_recursive_property_update": false,
            "thickness_projection_settings": {}
        }""")

        self.parameters.ValidateAndAssignDefaults(default_settings)

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts were provided for ShellThicknessControl. [ control name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", self.parameters["controlled_model_part_names"].GetStringArray(), False)
        self.model_part: 'Optional[Kratos.ModelPart]' = None

        self.filter = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), Kratos.THICKNESS, Kratos.Globals.DataLocation.Element, self.parameters["filter_settings"])

        self.output_all_fields = self.parameters["output_all_fields"].GetBool()
        self.physical_thicknesses = self.parameters["physical_thicknesses"].GetVector()

        self.thickness_projection = CreateProjection(parameters["thickness_projection_settings"], self.optimization_problem)

        self.consider_recursive_property_update = parameters["consider_recursive_property_update"].GetBool()

        num_phys_thick = len(self.physical_thicknesses)
        if num_phys_thick == 0:
            raise RuntimeError("The physical_thicknesses can not be empty, at least min and max should be provided.")
        self.filtered_thicknesses = [i for i, _ in enumerate(self.physical_thicknesses)]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        if not KratosOA.OptAppModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements, self.consider_recursive_property_update)
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        # initialize the filter
        self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
        self.filter.Initialize()

        # initialize the projections
        self.thickness_projection.SetProjectionSpaces(self.filtered_thicknesses, self.physical_thicknesses)

        physical_thickness_field = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(physical_thickness_field, Kratos.THICKNESS)

        # project backward the uniform physical control field and assign it to the control field
        self.physical_phi_field = self.thickness_projection.ProjectBackward(physical_thickness_field)

        # take the physical and control field the same
        self.control_phi_field = self.filter.UnfilterField(self.physical_phi_field)

        # now update the physical thickness field
        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        self.filter.Check()

    def Finalize(self) -> None:
        self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.THICKNESS]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_phi_field

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
            filtered_gradient = self.filter.BackwardFilterIntegratedField(projected_gradient)

            return filtered_gradient

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")

        update = new_control_field - self.control_phi_field
        if Kratos.Expression.Utils.NormL2(update) > 1e-15:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}.",False):
                # update the control thickness field
                self.control_phi_field = new_control_field
                # now update the physical field
                self._UpdateAndOutputFields(update)

                self.filter.Update()

                self.thickness_projection.Update()
                return True

        self.thickness_projection.Update()
        return False

    def _UpdateAndOutputFields(self, update: ContainerExpressionTypes) -> None:
        # filter the control field
        filtered_thickness_field_update = self.filter.ForwardFilterField(update)
        self.physical_phi_field = Kratos.Expression.Utils.Collapse(self.physical_phi_field + filtered_thickness_field_update)

        # project forward the filtered thickness field
        projected_filtered_thickness_field = self.thickness_projection.ProjectForward(self.physical_phi_field)
        # now update physical field
        KratosOA.PropertiesVariableExpressionIO.Write(projected_filtered_thickness_field, Kratos.THICKNESS)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(projected_filtered_thickness_field.GetContainer(), Kratos.THICKNESS)

        # compute and store projection derivatives for consistent filtering of the sensitivities
        self.projection_derivative_field = self.thickness_projection.ForwardProjectionGradient(self.physical_phi_field)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("physical_thickness", projected_filtered_thickness_field.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue("filtered_thickness_update", filtered_thickness_field_update.Clone(), overwrite=True)
            un_buffered_data.SetValue("control_thickness_phi", self.control_phi_field.Clone(), overwrite=True)
            un_buffered_data.SetValue("physical_thickness_phi", self.physical_phi_field.Clone(), overwrite=True)
            un_buffered_data.SetValue("projection_derivative", self.projection_derivative_field.Clone(), overwrite=True)

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {Kratos.THICKNESS.Name()}"
