import math
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import IntervalBounder
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"DataValuesControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DataValuesControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return DataValuesControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class DataValuesControl(Control):
    """Data values control

    This is a generic data values control which creates a control
    for the specified control variable. This does not do any filtering.

    """
    class __ElementValueControlHelper:
        clamper_type = KratosSI.ElementSmoothClamper
        data_location = Kratos.Globals.DataLocation.Element
        expression_type = Kratos.Expression.ElementExpression

        def ReadExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Read(exp, var)

        def WriteExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Write(exp, var)

    class __ConditionValueControlHelper:
        clamper_type = KratosSI.ConditionSmoothClamper
        data_location = Kratos.Globals.DataLocation.Condition
        expression_type = Kratos.Expression.ConditionExpression

        def ReadExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Read(exp, var)

        def WriteExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Write(exp, var)

    class __NodalValueControlHelper:
        clamper_type = KratosSI.NodeSmoothClamper
        data_location = Kratos.Globals.DataLocation.NodeNonHistorical
        expression_type = Kratos.Expression.NodalExpression

        def ReadExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Read(exp, var, False)

        def WriteExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Write(exp, var, False)

    class __NodalHistoricalValueControlHelper:
        clamper_type = KratosSI.NodeSmoothClamper
        data_location = Kratos.Globals.DataLocation.NodeHistorical
        expression_type = Kratos.Expression.NodalExpression

        def ReadExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Read(exp, var, True)

        def WriteExpression(self, exp, var):
            Kratos.Expression.VariableExpressionIO.Write(exp, var, True)

    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "control_variable_name"             : "",
            "control_variable_bounds"           : [0.0, 0.0],
            "container_type"                    : "",
            "output_all_fields"                 : false,
            "filter_settings"                   : {},
            "model_part_names": [
                {
                    "primal_model_part_name" : "PLEASE_PROVIDE_MODEL_PART_NAME",
                    "adjoint_model_part_name": "PLEASE_PROVIDE_MODEL_PART_NAME"
                }
            ]
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        allowed_container_types = {
            "nodal_historical": self.__NodalHistoricalValueControlHelper,
            "nodal_nonhistorical": self.__NodalValueControlHelper,
            "condition": self.__ConditionValueControlHelper,
            "element": self.__ElementValueControlHelper
        }
        container_type = parameters["container_type"].GetString()
        if container_type not in allowed_container_types.keys():
            raise RuntimeError(f"Unsupported container type = \"{container_type}\". Allowed container types are:\n\t" + "\n\t".join(list(allowed_container_types.keys())))
        else:
            self.helper = allowed_container_types[container_type]
            self.container_type_helper = self.helper()

        self.model : Kratos.Model = model
        self.parameters = parameters
        self.optimization_problem = optimization_problem

        self.output_all_fields = parameters["output_all_fields"].GetBool()
        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)
        if control_variable_type != "Double":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports double variables")
        self.controlled_physical_variable: SupportedSensitivityFieldVariableTypes = Kratos.KratosGlobals.GetVariable(control_variable_name)

        controlled_model_part_names: 'list[Kratos.Parameters]' = parameters["model_part_names"].values()
        if len(controlled_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for DataValuesControl. [ control name = \"{self.GetName()}\"]")

        self.primal_model_part_operation = ModelPartOperation(
                                                self.model,
                                                ModelPartOperation.OperationType.UNION,
                                                f"control_primal_{self.GetName()}",
                                                [param["primal_model_part_name"].GetString() for param in controlled_model_part_names],
                                                False)
        self.adjoint_model_part_operation = ModelPartOperation(
                                                self.model,
                                                ModelPartOperation.OperationType.UNION,
                                                f"control_adjoint_{self.GetName()}",
                                                [param["adjoint_model_part_name"].GetString() for param in controlled_model_part_names],
                                                False)

        # filter needs to be based on the primal model part
        # because, filter may keep pointers for the elements to get their center positions
        # for filtering. The adjoint model part may re-assign adjoint elements based on
        # the primal model part rendering the filter element pointers useless and to segfault.
        # hence filter is using the primal model part to get the locations.
        self.filter = FilterFactory(self.model, self.primal_model_part_operation.GetModelPartFullName(), self.controlled_physical_variable, self.container_type_helper.data_location, self.parameters["filter_settings"])

        self.primal_model_part: Optional[Kratos.ModelPart] = None
        self.adjoint_model_part: Optional[Kratos.ModelPart] = None

        control_variable_bounds = parameters["control_variable_bounds"].GetVector()

        # use the clamper in the unit interval
        self.interval_bounder = IntervalBounder(control_variable_bounds)
        self.clamper = self.container_type_helper.clamper_type(0, 1)

    def Initialize(self) -> None:
        self.primal_model_part = self.primal_model_part_operation.GetModelPart()
        self.adjoint_model_part = self.adjoint_model_part_operation.GetModelPart()

        # initialize the filter
        self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
        self.filter.Initialize()

        physical_field = self.GetPhysicalField()

        # get the phi field which is in [0, 1] range
        self.physical_phi_field = self.clamper.ProjectBackward(self.interval_bounder.GetBoundedExpression(physical_field))

        # compute the control phi field
        self.control_phi_field = self.filter.UnfilterField(self.physical_phi_field)

        self.physical_phi_derivative_field = self.clamper.CalculateForwardProjectionGradient(self.physical_phi_field) * self.interval_bounder.GetBoundGap()

        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        self.filter.Check()

    def Finalize(self) -> None:
        self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.controlled_physical_variable]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = self.container_type_helper.expression_type(self.primal_model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_phi_field

    def GetPhysicalField(self) -> ContainerExpressionTypes:
        physical_thickness_field = self.container_type_helper.expression_type(self.primal_model_part)
        self.container_type_helper.ReadExpression(physical_thickness_field, self.controlled_physical_variable)
        return physical_thickness_field

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        with TimeLogger("DataValuesControl::MapGradient", None, "Finished",False):
            keys = physical_gradient_variable_container_expression_map.keys()
            if len(keys) != 1:
                raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
            if self.controlled_physical_variable not in keys:
                raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {self.controlled_physical_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

            physical_gradient = physical_gradient_variable_container_expression_map[self.controlled_physical_variable]
            if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
                raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.primal_model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

            # dj/dE -> physical_gradient
            # dj/dphi = dj/dphysical * dphysical/dphi
            return self.filter.BackwardFilterIntegratedField(physical_gradient * self.physical_phi_derivative_field)

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.primal_model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        update = new_control_field - self.control_phi_field
        if not math.isclose(Kratos.Expression.Utils.NormL2(update), 0.0, abs_tol=1e-16):
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}.",False):
                # update the control thickness field
                self.control_phi_field = new_control_field
                # now update the physical field
                self._UpdateAndOutputFields(update)

                self.filter.Update()

                return True
        return False

    def _UpdateAndOutputFields(self, update: ContainerExpressionTypes) -> None:
        # filter the control field
        filtered_phi_field_update = self.filter.ForwardFilterField(update)
        self.physical_phi_field = Kratos.Expression.Utils.Collapse(self.physical_phi_field + filtered_phi_field_update)

        # project forward the filtered thickness field to get clamped physical field
        physical_field = self.interval_bounder.GetUnboundedExpression(self.clamper.ProjectForward(self.physical_phi_field))

        # now update physical field
        self.container_type_helper.WriteExpression(physical_field, self.controlled_physical_variable)

        # compute and store projection derivatives for consistent filtering of the sensitivities
        # this is dphi/dphysical -> physical_phi_derivative_field
        self.physical_phi_derivative_field = self.clamper.CalculateForwardProjectionGradient(self.physical_phi_field) * self.interval_bounder.GetBoundGap()

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue(f"{self.controlled_physical_variable.Name()}", physical_field.Clone(), overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue(f"{self.controlled_physical_variable.Name()}_physical_phi_update", filtered_phi_field_update.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"{self.controlled_physical_variable.Name()}_control_phi", self.control_phi_field.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"{self.controlled_physical_variable.Name()}_physical_phi", self.physical_phi_field.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"{self.controlled_physical_variable.Name()}_physical_phi_derivative", self.physical_phi_derivative_field.Clone(), overwrite=True)

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.adjoint_model_part_operation.GetModelPartFullName()}, control variable = {self.controlled_physical_variable.Name()}"