from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.filtering.filtering_factory import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"ShapeControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"ShapeControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return ShapeControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class ShapeControl(Control):
    """Node-based shape control using implicit and explicit Vertex Morphing techniques

    This is filtering-based discrete shape control which parametrizes and controls
    discrete shell and solid geometries.

    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [],
            "filter_settings"            : {},
            "use_mesh_motion_solver"     : false,
            "mesh_motion_settings"       : {},
            "output_all_fields"          : false
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_all_fields = parameters["output_all_fields"].GetBool()

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", parameters["controlled_model_part_names"].GetStringArray(), False)
        self.model_part: 'Optional[Kratos.ModelPart]' = None

        self.filter = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), KratosOA.SHAPE, Kratos.Globals.DataLocation.NodeHistorical, parameters["filter_settings"])

        # set the mesh motion solver if required
        self.use_mesh_motion_solver = parameters["use_mesh_motion_solver"].GetBool()
        if self.use_mesh_motion_solver:
            from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
            self.mesh_motion_solver = MeshMovingAnalysis(self.model, parameters["mesh_motion_settings"])

    @time_decorator(methodName="Initialize")
    def Initialize(self) -> None:
        self.filter.Initialize()
        if self.use_mesh_motion_solver:
            self.mesh_motion_solver.Initialize()

    def Check(self) -> None:
        self.filter.Check()

    def Finalize(self) -> None:
        self.filter.Finalize()
        if self.use_mesh_motion_solver:
            self.mesh_motion_solver.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, [0,0,0])
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field

    def GetPhysicalField(self) -> ContainerExpressionTypes:
        physical_shape_field = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(physical_shape_field, Kratos.Configuration.Initial)
        return physical_shape_field

    @time_decorator(methodName="MapGradient")
    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if KratosOA.SHAPE not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {KratosOA.SHAPE.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        physical_gradient = physical_gradient_variable_container_expression_map[KratosOA.SHAPE]
        if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
            raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.filter_model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

        # now filter the field
        filtered_gradient = self.filter.FilterIntegratedField(physical_gradient)

        return filtered_gradient

    @time_decorator(methodName="Update")
    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.filter_model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")
        if Kratos.Expression.Utils.NormL2(self.control_field - new_control_field) > 1e-15:
            # update the control SHAPE field
            control_update = new_control_field - self.control_field
            self.control_field = new_control_field
            # now update the physical field
            self._UpdateAndOutputFields(control_update)
            # now update the filter
            self.filter.Update()
            return True
        return False

    def _UpdateAndOutputFields(self, control_update: ContainerExpressionTypes) -> None:
        # compute the shape update
        shape_update = self.filter.FilterField(control_update)

        self._UpdateMesh(shape_update)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("shape_update", shape_update.Clone(), overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue("shape_control", self.control_field.Clone(), overwrite=True)
            un_buffered_data.SetValue("shape_control_update", control_update.Clone(), overwrite=True)

    def _UpdateMesh(self, shape_update: ContainerExpressionTypes) -> None:
        if self.use_mesh_motion_solver:
            mm_computing_model_part: Kratos.ModelPart = self.mesh_motion_solver._GetSolver().GetComputingModelPart()
            time_before_update = mm_computing_model_part.ProcessInfo.GetValue(Kratos.TIME)
            step_before_update = mm_computing_model_part.ProcessInfo.GetValue(Kratos.STEP)
            delta_time_before_update = mm_computing_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)

            # Reset step/time iterators such that they match the current iteration after calling RunSolutionLoop (which internally calls CloneTimeStep)
            mm_computing_model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_update-1)
            mm_computing_model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_update-1)
            mm_computing_model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0)

            # assign the boundary_displacements
            Kratos.Expression.VariableExpressionIO.Write(shape_update, Kratos.MESH_DISPLACEMENT, True)

            # solve for the volume mesh displacements
            if not self.mesh_motion_solver.time < self.mesh_motion_solver.end_time:
                self.mesh_motion_solver.end_time += 1
            self.mesh_motion_solver.RunSolutionLoop()

            # update the coords
            mesh_displacement_field = Kratos.Expression.NodalExpression(mm_computing_model_part)
            physical_shape_field = Kratos.Expression.NodalExpression(mm_computing_model_part)
            Kratos.Expression.VariableExpressionIO.Read(mesh_displacement_field, Kratos.MESH_DISPLACEMENT, True)
            Kratos.Expression.NodalPositionExpressionIO.Read(physical_shape_field, Kratos.Configuration.Initial)

            Kratos.Expression.NodalPositionExpressionIO.Write(mesh_displacement_field + physical_shape_field, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(mesh_displacement_field + physical_shape_field, Kratos.Configuration.Current)

            mm_computing_model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_update)
            mm_computing_model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_update)
            mm_computing_model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, delta_time_before_update)
        else:
            physical_field = self.GetPhysicalField()
            Kratos.Expression.NodalPositionExpressionIO.Write(shape_update + physical_field, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(shape_update + physical_field, Kratos.Configuration.Current)

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, filter = {self.filter}, control variable = {KratosOA.SHAPE.Name()}"
