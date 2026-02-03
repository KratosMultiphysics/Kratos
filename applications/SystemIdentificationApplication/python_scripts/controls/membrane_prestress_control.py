import math
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionBoundingManager
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.SystemIdentificationApplication.controls.material_properties_control import MaterialPropertiesControl

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"MembranePrestressControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MembranePrestressControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return MembranePrestressControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class MembranePrestressControl(MaterialPropertiesControl):
    """Material properties control

    This is a membrane presress in material properties which creates a control
    for the specified control variable. This does not do any filtering.

    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super(MaterialPropertiesControl, self).__init__(name)

        default_settings = Kratos.Parameters("""{
            "control_variable_name"             : "",
            "control_variable_bounds"           : [0.0, 0.0],
            "output_all_fields"                 : false,
            "consider_recursive_property_update": false,
            "filter_settings"                   : {},
            "model_part_names": [
                {
                    "primal_model_part_name" : "PLEASE_PROVIDE_MODEL_PART_NAME",
                    "adjoint_model_part_name": "PLEASE_PROVIDE_MODEL_PART_NAME"
                }
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.parameters = parameters
        self.optimization_problem = optimization_problem

        self.output_all_fields = parameters["output_all_fields"].GetBool()
        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)
        if control_variable_type != "Array":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports Array1DVariable3 variables")
        self.controlled_physical_variable: SupportedSensitivityFieldVariableTypes = Kratos.KratosGlobals.GetVariable(control_variable_name)

        controlled_model_part_names: 'list[Kratos.Parameters]' = parameters["model_part_names"].values()
        if len(controlled_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for MembranePrestressControl. [ control name = \"{self.GetName()}\"]")

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

        self.consider_recursive_property_update = parameters["consider_recursive_property_update"].GetBool()


        # filter needs to be based on the primal model part
        # because, filter may keep pointers for the elements to get their center positions
        # for filtering. The adjoint model part may re-assign adjoint elements based on
        # the primal model part rendering the filter element pointers useless and to segfault.
        # hence filter is using the primal model part to get the locations.
        self.filter = FilterFactory(self.model, self.primal_model_part_operation.GetModelPartFullName(), self.controlled_physical_variable, Kratos.Globals.DataLocation.Element, self.parameters["filter_settings"])
        #print(self.filter)
        self.primal_model_part: Optional[Kratos.ModelPart] = None
        self.adjoint_model_part: Optional[Kratos.ModelPart] = None

        control_variable_bounds = parameters["control_variable_bounds"].GetVector()

        # use the clamper in the unit interval
        self.interval_bounder = ExpressionBoundingManager(control_variable_bounds)
        self.clamper = KratosSI.ElementSmoothClamper(0, 1)

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.primal_model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, Kratos.Array3(0.0))
        return field

    def GetPhysicalField(self) -> ContainerExpressionTypes:
        physical_thickness_field = Kratos.Expression.ElementExpression(self.primal_model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(physical_thickness_field, KratosSM.PRESTRESS_VECTOR)
        return physical_thickness_field

    def _UpdateAndOutputFields(self, update: ContainerExpressionTypes) -> None:
        #print("In update")
        #print("update ", update.Evaluate())
        # filter the control field
        filtered_phi_field_update = self.filter.ForwardFilterField(update)
        #print("filtered_phi_field_update ", filtered_phi_field_update.Evaluate())
        self.physical_phi_field = Kratos.Expression.Utils.Collapse(self.physical_phi_field + filtered_phi_field_update)
        #print("self.physical_phi_field ", self.physical_phi_field.Evaluate())

        # project forward the filtered thickness field to get clamped physical field
        physical_field = self.interval_bounder.GetUnboundedExpression(self.clamper.ProjectForward(self.physical_phi_field))
        #print(self.controlled_physical_variable.Name(), " physical_field ", physical_field.Evaluate())
        # now update physical field
        KratosOA.PropertiesVariableExpressionIO.Write(physical_field, KratosSM.PRESTRESS_VECTOR)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(physical_field.GetContainer(), KratosSM.PRESTRESS_VECTOR)

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

    def _GetControlPrefix(self) -> str:
        return "MembranePrestressControl"
    