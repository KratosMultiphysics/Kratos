from typing import Optional
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.filtering.filter import Filter
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"ExternalResponseFunctionControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"ExternalResponseFunctionControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return ExternalResponseFunctionControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class ExternalResponseFunctionControl(Control):
    def __init__(self, control_name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(control_name)

        default_parameters = Kratos.Parameters("""{
            "controlled_model_part_names": [],
            "design_variable"            : "CUSTOM_DESIGN_VARIABLE",
            "container_type"             : "nodal_nonhistorical",
            "use_filtering"              : false,
            "filter_settings"            : {},
            "output_all_fields"          : false
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.optimization_problem = optimization_problem

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", parameters["controlled_model_part_names"].GetStringArray(), False)
        self.model_part: 'Optional[Kratos.ModelPart]' = None

        self.design_variable: SupportedSensitivityFieldVariableTypes = Kratos.KratosGlobals.GetVariable(parameters["design_variable"].GetString())
        self.output_all_fields = parameters["output_all_fields"].GetBool()

        self.container_type = parameters["container_type"].GetString()

        if parameters["use_filtering"].GetBool():
            self.filter: 'Optional[Filter]' = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), self.design_variable, self.__GetDataLocation(), parameters["filter_settings"])
        else:
            self.filter = None


    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        ComponentDataView(self, self.optimization_problem).SetDataBuffer(1)

        # initialize the filter
        if self.filter:
            self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
            self.filter.Initialize()

        self.physical_phi_field = self.__CreateExpression(self.model_part)
        self.__ReadExpression(self.physical_phi_field, self.design_variable)

        if self.filter:
            # take the physical and control field the same
            self.control_phi_field = self.filter.UnfilterField(self.physical_phi_field)
        else:
            self.control_phi_field = self.physical_phi_field.Clone()

        # now update the physical thickness field
        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        if self.filter:
            self.filter.Check()

    def Finalize(self) -> None:
        if self.filter:
            self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.design_variable]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        exp = self.__CreateExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetDataToZero(exp, self.design_variable)
        return exp

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_phi_field

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        with TimeLogger("ExternalResponseFunctionControl::MapGradient", None, "Finished",False):
            keys = physical_gradient_variable_container_expression_map.keys()
            if len(keys) != 1:
                raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
            if self.design_variable not in keys:
                raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {self.design_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

            physical_gradient = physical_gradient_variable_container_expression_map[self.design_variable]
            if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
                raise RuntimeError(f"Gradients for the required container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

            if self.filter:
                # now filter the field
                filtered_gradient = self.filter.BackwardFilterIntegratedField(physical_gradient)
            else:
                filtered_gradient = physical_gradient.Clone()

            return filtered_gradient

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")

        update = new_control_field - self.control_phi_field
        if Kratos.Expression.Utils.NormL2(update) > 1e-15:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}.",False):
                # update the control thickness field
                self.control_phi_field = new_control_field
                # now update the physical field
                self._UpdateAndOutputFields(update)

                if self.filter:
                    self.filter.Update()
                return True
        return False

    def _UpdateAndOutputFields(self, update: ContainerExpressionTypes) -> None:
        if self.filter:
            # filter the control field
            filtered_update = self.filter.ForwardFilterField(update)
        else:
            filtered_update = update.Clone()

        self.physical_phi_field = Kratos.Expression.Utils.Collapse(self.physical_phi_field + filtered_update)
        self.__WriteExpression(self.physical_phi_field, self.design_variable)

        # add the values of the control to the optimization problem.
        component_data_view = ComponentDataView(self, self.optimization_problem)
        for i, v in enumerate(self.physical_phi_field.Evaluate()):
            component_data_view.GetBufferedData().SetValue(f"{self.model_part.FullName()}_point_{i+1}", float(v))

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}", self.physical_phi_field.Clone(), overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}_update", update.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}_filtered_update", filtered_update.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}_control_phi", self.control_phi_field.Clone(), overwrite=True)

    def __CreateExpression(self, model_part: Kratos.ModelPart) -> ContainerExpressionTypes:
        exp_dict = {
            "nodal_historical": lambda : Kratos.Expression.NodalExpression(model_part),
            "nodal_nonhistorical": lambda : Kratos.Expression.NodalExpression(model_part),
            "condition": lambda : Kratos.Expression.NodalExpression(model_part),
            "element": lambda : Kratos.Expression.NodalExpression(model_part),
        }
        return exp_dict[self.container_type]()

    def __ReadExpression(self, expression: ContainerExpressionTypes, variable: SupportedSensitivityFieldVariableTypes) -> None:
        exp_dict = {
            "nodal_historical": lambda : Kratos.Expression.VariableExpressionIO.Read(expression, variable, True),
            "nodal_nonhistorical": lambda : Kratos.Expression.VariableExpressionIO.Read(expression, variable, False),
            "condition": lambda : Kratos.Expression.VariableExpressionIO.Read(expression, variable),
            "element": lambda :Kratos.Expression.VariableExpressionIO.Read(expression, variable),
        }
        return exp_dict[self.container_type]()

    def __WriteExpression(self, expression: ContainerExpressionTypes, variable: SupportedSensitivityFieldVariableTypes) -> None:
        exp_dict = {
            "nodal_historical": lambda : Kratos.Expression.VariableExpressionIO.Write(expression, variable, True),
            "nodal_nonhistorical": lambda : Kratos.Expression.VariableExpressionIO.Write(expression, variable, False),
            "condition": lambda : Kratos.Expression.VariableExpressionIO.Write(expression, variable),
            "element": lambda :Kratos.Expression.VariableExpressionIO.Write(expression, variable),
        }
        return exp_dict[self.container_type]()

    def __GetDataLocation(self) -> Kratos.Globals.DataLocation:
        exp_dict = {
            "nodal_historical": Kratos.Globals.DataLocation.NodeHistorical,
            "nodal_nonhistorical": Kratos.Globals.DataLocation.NodeNonHistorical,
            "condition": Kratos.Globals.DataLocation.Condition,
            "element": Kratos.Globals.DataLocation.Element
        }
        return exp_dict[self.container_type]


