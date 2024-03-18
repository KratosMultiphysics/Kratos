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
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"VertexMorphingShapeControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"VertexMorphingShapeControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return VertexMorphingShapeControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class VertexMorphingShapeControl(Control):
    """Node-based shape control using implicit and explicit Vertex Morphing techniques

    This is filtering-based discrete shape control which parametrizes and controls
    discrete shell and solid geometries.

    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [],
            "filter_settings"            : {},
            "mesh_motion_solver_type"    : "none",
            "skin_model_part_names"      : [],
            "mesh_motion_solver_settings": {},
            "output_all_fields"          : false
        }""")

        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.output_all_fields = self.parameters["output_all_fields"].GetBool()

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for VertexMorphingShapeControl. [ control name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'Optional[Kratos.ModelPart]' = None

        self.filter = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), KratosOA.SHAPE, Kratos.Globals.DataLocation.NodeHistorical, self.parameters["filter_settings"])

        # now setup mesh motion
        self.__mesh_moving_analysis: 'Optional[MeshMovingAnalysis]' = None
        self.__boundary_model_part: 'Optional[Kratos.ModelPart]' = None

        self.mesh_motion_solver_type = self.parameters["mesh_motion_solver_type"].GetString()
        supported_mesh_motion_solver_types = ["mesh_moving_analysis", "filter_based", "none"]
        if self.mesh_motion_solver_type not in supported_mesh_motion_solver_types:
            raise RuntimeError(f"Unsupported mesh_motion_solver_type = \"{self.mesh_motion_solver_type}\" requested. Supported types are:\n\t" + "\n\t".join(supported_mesh_motion_solver_types))

        if self.mesh_motion_solver_type == "mesh_moving_analysis":
            # import is done here, so that, if it is not required to use mesh_motion,
            # then it is not loaded, hence not even required to have MeshMovingApplication
            # compiled.
            from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
            self.__mesh_moving_analysis = MeshMovingAnalysis(self.model, self.parameters["mesh_motion_solver_settings"])

    @time_decorator(methodName="GetName")
    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        # initialize filters
        self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
        self.filter.Initialize()

        self.control_field = self.filter.UnfilterField(self.GetPhysicalField())

        # initialize_mesh_motion if defined.
        if not self.__mesh_moving_analysis is None:
            self.__mesh_moving_analysis.Initialize()

            def recursively_increase_usage(model_part: Kratos.ModelPart) -> None:
                if not model_part.Has(KratosOA.NUMBER_OF_SOLVERS_USING_NODES):
                    model_part[KratosOA.NUMBER_OF_SOLVERS_USING_NODES] = 0
                model_part[KratosOA.NUMBER_OF_SOLVERS_USING_NODES] += 1

                for sub_model_part in model_part.SubModelParts:
                    recursively_increase_usage(sub_model_part)

            # we increase the usage of all the sub-model parts because, mesh moving works with
            # root model parts, and all the sub model parts are suing the root model parts. Hence
            # mesh moving will affect all the shared sub-model parts
            recursively_increase_usage(self.__mesh_moving_analysis._GetSolver().GetComputingModelPart())

        # now update
        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        self.filter.Check()

    def Finalize(self) -> None:
        self.filter.Finalize()
        if not self.__mesh_moving_analysis is None:
            self.__mesh_moving_analysis.Finalize()
            del self.__mesh_moving_analysis

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.NodalExpression(self.__GetPhysicalModelPart())
        Kratos.Expression.LiteralExpressionIO.SetData(field, [0,0,0])
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field

    def GetPhysicalField(self) -> ContainerExpressionTypes:
        physical_shape_field = Kratos.Expression.NodalExpression(self.__GetPhysicalModelPart())
        Kratos.Expression.NodalPositionExpressionIO.Read(physical_shape_field, Kratos.Configuration.Initial)
        return physical_shape_field

    @time_decorator(methodName="GetName")
    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if KratosOA.SHAPE not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {KratosOA.SHAPE.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        physical_gradient = physical_gradient_variable_container_expression_map[KratosOA.SHAPE]
        if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
            raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

        filtered_gradient = self.filter.BackwardFilterIntegratedField(KratosOA.ExpressionUtils.ExtractData(physical_gradient, self.model_part))

        return KratosOA.ExpressionUtils.ExtractData(filtered_gradient, self.__GetPhysicalModelPart())

    @time_decorator(methodName="GetName")
    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")
        if Kratos.Expression.Utils.NormL2(self.control_field - new_control_field) > 1e-15:
            # update the control SHAPE field
            control_update = new_control_field - self.control_field
            self.control_field = new_control_field
            # now update the physical field
            self._UpdateAndOutputFields(control_update)

            self.filter.Update()
            return True
        return False

    def _UpdateAndOutputFields(self, control_update: ContainerExpressionTypes) -> None:
        # compute the shape update
        if self.mesh_motion_solver_type == "filter_based":
            shape_update = self.filter.ForwardFilterField(control_update)
        else:
            shape_update = self.filter.ForwardFilterField(KratosOA.ExpressionUtils.ExtractData(control_update, self.model_part))

        # now update the shape
        self._UpdateMesh(shape_update)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("shape_update", shape_update.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue("shape_control", self.control_field.Clone(),overwrite=True)
            un_buffered_data.SetValue("shape_control_update", control_update.Clone(),overwrite=True)

    def _UpdateMesh(self, shape_update: ContainerExpressionTypes) -> None:
        if self.mesh_motion_solver_type != "none":
            mesh_displacement_field = Kratos.Expression.NodalExpression(self.model_part_operation.GetRootModelPart())

            if not self.__mesh_moving_analysis is None:
                mm_computing_model_part: Kratos.ModelPart = self.__mesh_moving_analysis._GetSolver().GetComputingModelPart()
                time_before_update = mm_computing_model_part.ProcessInfo.GetValue(Kratos.TIME)
                step_before_update = mm_computing_model_part.ProcessInfo.GetValue(Kratos.STEP)
                delta_time_before_update = mm_computing_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)

                # Reset step/time iterators such that they match the current iteration after calling RunSolutionLoop (which internally calls CloneTimeStep)
                mm_computing_model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_update-1)
                mm_computing_model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_update-1)
                mm_computing_model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0)

                # first reset the MESH_DISPLACEMENTS
                Kratos.VariableUtils().SetHistoricalVariableToZero(Kratos.MESH_DISPLACEMENT, mm_computing_model_part.Nodes)

                # assign the design surface MESH_DISPLACEMENTS
                Kratos.Expression.VariableExpressionIO.Write(shape_update, Kratos.MESH_DISPLACEMENT, True)

                # since we know that the nodes in the shape_update model part needs to be fixed to correctly
                # compute the MESH_DISPLACEMENT. We Only do that in here by first freeing all the MESH_DISPLACEMENT
                # dofs, then followed by fixing the shape_update model part nodes' MESH_DISPLACEMENTs. If further
                # boundaries needs fixing, then they need to be specified in the mesh_motion_solver settings
                # using fix_vector_variable_process.
                mesh_displacement_comps = [Kratos.MESH_DISPLACEMENT_X, Kratos.MESH_DISPLACEMENT_Y, Kratos.MESH_DISPLACEMENT_Z]
                for i_comp, mesh_displacement_var in enumerate(mesh_displacement_comps):
                    # first free the mesh displacement component
                    Kratos.VariableUtils().ApplyFixity(mesh_displacement_var, False, mm_computing_model_part.Nodes)
                    # now fix the design surface component
                    Kratos.VariableUtils().ApplyFixity(mesh_displacement_var, True, shape_update.GetModelPart().Nodes)
                    # now fix boundary condition components
                    for model_part in self.filter.GetBoundaryConditions()[i_comp]:
                        Kratos.VariableUtils().ApplyFixity(mesh_displacement_var, True, model_part.Nodes)

                # solve for the volume mesh displacements
                if not self.__mesh_moving_analysis.time < self.__mesh_moving_analysis.end_time:
                    self.__mesh_moving_analysis.end_time += 1
                self.__mesh_moving_analysis.RunSolutionLoop()

                Kratos.Expression.VariableExpressionIO.Read(mesh_displacement_field, Kratos.MESH_DISPLACEMENT, True)

                mm_computing_model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_update)
                mm_computing_model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_update)
                mm_computing_model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, delta_time_before_update)
            else:
                # filter based mesh motion is used. Hence read the value from the filter. then the shape update is the mesh displacement
                mesh_displacement_field = shape_update.Clone()

            # a mesh motion is solved (either by mesh motion solver or the filter), and MESH_DISPLACEMENT variable is filled
            # Hence, update the coords
            physical_field = Kratos.Expression.NodalExpression(self.model_part_operation.GetRootModelPart())

            Kratos.Expression.NodalPositionExpressionIO.Read(physical_field, Kratos.Configuration.Initial)

            Kratos.Expression.NodalPositionExpressionIO.Write(mesh_displacement_field + physical_field, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(mesh_displacement_field + physical_field, Kratos.Configuration.Current)
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Applied MESH_DISPLACEMENT to the mesh.")
        else:
            # no mesh motion is preferred, hence only updating the design field
            physical_field = self.GetPhysicalField()
            Kratos.Expression.NodalPositionExpressionIO.Write(shape_update + physical_field, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(shape_update + physical_field, Kratos.Configuration.Current)

    def __GetPhysicalModelPart(self) -> Kratos.ModelPart:
        if self.mesh_motion_solver_type == "none":
            # this control does not require mesh_motion, hence this does not
            # have any influence on the volume domain. So the physical model part
            # is the controlled model part
            return self.model_part
        else:
            # this is the case where mesh motion is required, that means this control will control the
            # whole domain by moving its internal nodes because, the controlled model part represents
            # a surface of the volume domain. hence the physical field should be the whole
            # model part
            return self.model_part.GetRootModelPart()

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = SHAPE ]"