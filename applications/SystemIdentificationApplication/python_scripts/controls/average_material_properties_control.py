import math
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
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
        raise RuntimeError(f"AverageMaterialPropertiesControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"AverageMaterialPropertiesControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return AverageMaterialPropertiesControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class AverageMaterialPropertiesControl(MaterialPropertiesControl):
    """Material properties control

    This is a generic material properties control which creates a control
    for the specified control variable. This does not do any filtering.

    """
    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        with TimeLogger("ShellThicknessControl::MapGradient", None, "Finished",False):
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
            return self.filter.BackwardFilterField(physical_gradient * self.physical_phi_derivative_field)
        
    def _GetControlPrefix(self) -> str:
        return "AverageMaterialPropertiesControl"