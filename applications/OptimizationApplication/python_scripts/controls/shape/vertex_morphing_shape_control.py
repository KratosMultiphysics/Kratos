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
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

try:
    from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
    mesh_motion_solver_available = True
except ImportError:
    mesh_motion_solver_available = False

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"VertexMorphingShapeControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"VertexMorphingShapeControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return VertexMorphingShapeControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

def GetImplicitFilterParameters(filter_model_part_name: str, parameters: Kratos.Parameters):
    filter_radius = parameters["filter_settings"]["radius"].GetDouble()
    filter_type = parameters["filter_settings"]["type"].GetString()
    linear_solver_settings = parameters["filter_settings"]["linear_solver_settings"]
    implicit_vector_filter_parameters = Kratos.Parameters("""
     {
        "solver_settings" : {
             "domain_size"          : 3,
             "echo_level"           : 0,
             "filter_type"          : "general_vector",
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

    if filter_type == "bulk_surface_implicit":
        implicit_vector_filter_parameters["solver_settings"]["filter_type"].SetString("bulk_surface_shape")

    implicit_vector_filter_parameters["solver_settings"].AddDouble("filter_radius",filter_radius)
    implicit_vector_filter_parameters["solver_settings"].AddString("model_part_name",filter_model_part_name)
    if linear_solver_settings.Has("solver_type"):
        implicit_vector_filter_parameters["solver_settings"].AddValue("linear_solver_settings", linear_solver_settings)

    for model_part_name, direction_params in parameters["fixed_model_parts"].items():
        auto_process_settings = Kratos.Parameters(
            """
            {
                "python_module" : "fix_vector_variable_process",
                "kratos_module" : "KratosMultiphysics",
                "help"          : "This process fixes the selected components of a given vector variable without modifying the value of the variable.",
                "process_name"  : "FixVectorVariableProcess",
                "Parameters"    : {
                    "model_part_name" : \""""+model_part_name+"""\",
                    "variable_name"   : "HELMHOLTZ_VECTOR",
                    "constrained"          : [false,false,false]
                }
            }
            """)
        auto_process_settings["Parameters"]["constrained"] = direction_params
        implicit_vector_filter_parameters["processes"]["boundary_conditions_process_list"].Append(auto_process_settings)

    return implicit_vector_filter_parameters

def GetMeshMotionParameters(parameters: Kratos.Parameters, main_model_part_name, moving_sub_model_part_name):
    default_mesh_motion_settings = Kratos.Parameters("""
    {
        "apply_mesh_solver" : true,
        "solver_settings" : {
            "domain_size"     : 3,
            "echo_level"      : 0,
            "solver_type"     : "structural_similarity",
            "model_part_name" : "NONE",
            "model_import_settings"              : {
                "input_type"     : "use_input_model_part"
            },
            "time_stepping" : {
                "time_step"       : 1.0
            },
            "compute_reactions"                : false,
            "calculate_mesh_velocity"          : false
        },
        "processes" : {
            "boundary_conditions_process_list" : []
        }
    }""")

    parameters.ValidateAndAssignDefaults(default_mesh_motion_settings)

    parameters["solver_settings"].AddString("model_part_name",main_model_part_name)

    auto_process_settings = Kratos.Parameters(
    """
    {
        "python_module" : "fix_vector_variable_process",
        "kratos_module" : "KratosMultiphysics",
        "help"          : "This process fixes the selected components of a given vector variable without modifying the value of the variable.",
        "process_name"  : "FixVectorVariableProcess",
        "Parameters"    : {
            "model_part_name"      : \""""+str(moving_sub_model_part_name)+"""\",
            "variable_name"        : "MESH_DISPLACEMENT",
            "constrained"          : [true,true,true]
        }
    }
    """)

    parameters["processes"]["boundary_conditions_process_list"].Append(auto_process_settings)

    problem_data = Kratos.Parameters("""{
        "echo_level"    : 0,
        "start_time"    : 0.0,
        "end_time"      : 1.0,
        "parallel_type" : "OpenMP"
    }""")

    parameters.AddValue("problem_data", problem_data)

    return parameters

class VertexMorphingShapeControl(Control):
    """Node-based shape control using implicit and explicit Vertex Morphing techniques

    This is filtering-based discrete shape control which parametrizes and controls
    discrete shell and solid geometeries.

    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [],
            "filter_settings": {
                "type" : "surface_explicit",
                "filter_function_type" : "linear",
                "damping_function_type" : "sigmoidal",
                "radius": 0.000000000001,
                "max_nodes_in_filter_radius": 1000,
                "linear_solver_settings" : {}
            },
            "mesh_motion_settings" : {},
            "output_all_fields": false,
            "fixed_model_parts": {},
            "utilities": []
        }""")

        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.parameters["filter_settings"].ValidateAndAssignDefaults(default_settings["filter_settings"])

        self.output_all_fields = self.parameters["output_all_fields"].GetBool()
        self.filter_type = self.parameters["filter_settings"]["type"].GetString()

        self.supported_filter_types = ["surface_implicit","surface_implicit_with_mesh_motion","bulk_surface_implicit","surface_explicit","surface_explicit_with_mesh_motion"]
        if not self.filter_type in self.supported_filter_types:
            raise RuntimeError(f"The specified filter type \"{self.filter_type}\" is not supported. Followings are the variables:\n\t" + "\n\t".join([k for k in self.supported_filter_types]))

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for VertexMorphingShapeControl. [ control name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'Optional[Kratos.ModelPart]' = None

        self.filter = None
        self.is_filter_implicit = False
        if "implicit" in self.filter_type:
            self.is_filter_implicit = True
            helmholtz_settings = GetImplicitFilterParameters(self.model_part_operation.GetModelPartFullName(), self.parameters)
            self.filter = HelmholtzAnalysis(self.model, helmholtz_settings)

        # now setup mesh motion
        self.SetupMeshMotion()

        self.is_initialized = False

    @time_decorator(methodName="GetName")
    def Initialize(self) -> None:

        if self.filter_type == "bulk_surface_implicit":
            self.filter_model_part = self.model_part_operation.GetModelPart().GetRootModelPart()
        else:
            self.filter_model_part = self.model_part_operation.GetModelPart()

        # initialize filters
        if not self.is_filter_implicit:
            # make sure that explicit filter model part only has conditions
            if len(self.filter_model_part.Elements) > 0:
                raise RuntimeError("VertexMorphingShapeControl with explicit filtering only allows model parts with conditions")
            # now damping and creating the filter
            fixed_model_parts_names = list(self.parameters["fixed_model_parts"].keys())
            if len(fixed_model_parts_names)>0:
                self.fixed_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", fixed_model_parts_names, False)
                self.fixed_model_part = self.fixed_model_part_operation.GetModelPart()
                self.filter = KratosOA.NodalExplicitFilter(self.filter_model_part, self.fixed_model_part, self.parameters["filter_settings"]["filter_function_type"].GetString(),
                                                           self.parameters["filter_settings"]["damping_function_type"].GetString(),
                                                           self.parameters["filter_settings"]["max_nodes_in_filter_radius"].GetInt())
            else:
                self.filter = KratosOA.NodalExplicitFilter(self.filter_model_part, self.parameters["filter_settings"]["filter_function_type"].GetString(),
                                                           self.parameters["filter_settings"]["max_nodes_in_filter_radius"].GetInt())

            filter_radius_field = Kratos.Expression.NodalExpression(self.filter_model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(filter_radius_field, self.parameters["filter_settings"]["radius"].GetDouble())
            self.filter.SetFilterRadius(filter_radius_field)
            self.control_field = Kratos.Expression.NodalExpression(self.filter_model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(self.control_field,[0,0,0])
        else:
            self.filter.Initialize()
            self._SetFixedModelPartValues(False)
            self.control_field = self.filter.UnFilterField(self.GetPhysicalField())
            self._SetFixedModelPartValues(True)

        # initialize_mesh_motion
        if self.has_mesh_motion:
            self._mesh_moving_analysis.Initialize()

        # now update
        self._UpdateAndOutputFields(self.GetEmptyField())
        self.is_initialized = True

    def SetupMeshMotion(self):
        if "mesh_motion" in self.filter_type:
            self.has_mesh_motion = True
            self.mesh_motion_params = GetMeshMotionParameters(self.parameters["mesh_motion_settings"],self.model_part_operation.GetRootModelPart().Name,self.model_part_operation.GetModelPartFullName())
            self._mesh_moving_analysis = MeshMovingAnalysis(self.model, self.mesh_motion_params)
        else:
            self.has_mesh_motion = False

    def Check(self) -> None:
        if self.is_filter_implicit:
            return self.filter.Check()

    def Finalize(self) -> None:
        if self.is_filter_implicit:
            self.filter.Finalize()
        if self.has_mesh_motion:
            self._mesh_moving_analysis.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.NodalExpression(self.filter_model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, [0,0,0])
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field

    def GetPhysicalField(self) -> ContainerExpressionTypes:
        physical_shape_field = Kratos.Expression.NodalExpression(self.filter_model_part)
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
            raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.filter_model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

        # now filter the field
        filtered_gradient = self.filter.FilterIntegratedField(physical_gradient)

        return filtered_gradient

    @time_decorator(methodName="GetName")
    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.filter_model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")
        if Kratos.Expression.Utils.NormL2(self.control_field - new_control_field) > 1e-15:
            # update the control SHAPE field
            control_update = new_control_field - self.control_field
            self.control_field = new_control_field
            # now update the physical field
            self._UpdateAndOutputFields(control_update)
            # update explicit filter
            if not self.is_filter_implicit:
                self.filter.Update()
            return True
        return False

    def _UpdateMesh(self, shape_update) -> None:
        if self.has_mesh_motion:
            mm_computing_model_part = self._mesh_moving_analysis._GetSolver().GetComputingModelPart()
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
            if not self._mesh_moving_analysis.time < self._mesh_moving_analysis.end_time:
                self._mesh_moving_analysis.end_time += 1
            self._mesh_moving_analysis.RunSolutionLoop()

            # update the coords
            mesh_disp_field = Kratos.Expression.NodalExpression(mm_computing_model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(mesh_disp_field, [0,0,0])
            Kratos.Expression.VariableExpressionIO.Read(mesh_disp_field, Kratos.MESH_DISPLACEMENT, True)
            physical_shape_field = Kratos.Expression.NodalExpression(mm_computing_model_part)
            Kratos.Expression.NodalPositionExpressionIO.Read(physical_shape_field, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(mesh_disp_field + physical_shape_field, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(mesh_disp_field + physical_shape_field, Kratos.Configuration.Current)

            mm_computing_model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_update)
            mm_computing_model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_update)
            mm_computing_model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, delta_time_before_update)
        else:
            Kratos.Expression.NodalPositionExpressionIO.Write(shape_update + self.GetPhysicalField(), Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(shape_update + self.GetPhysicalField(), Kratos.Configuration.Current)

    def _UpdateAndOutputFields(self, control_update) -> None:
        # compute the shape update
        shape_update = self.filter.FilterField(control_update)

        # now update the shape
        self._UpdateMesh(shape_update)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("shape_update", shape_update.Clone(),overwrite=self.is_initialized)
        if self.output_all_fields:
            un_buffered_data.SetValue("shape_control", self.control_field.Clone(),overwrite=self.is_initialized)
            un_buffered_data.SetValue("shape_control_update", control_update.Clone(),overwrite=self.is_initialized)

    def _SetFixedModelPartValues(self, is_backward = False) -> None:
        for model_part_name in self.parameters["fixed_model_parts"].keys():
            field =  Kratos.Expression.NodalExpression(self.model[model_part_name])
            Kratos.Expression.LiteralExpressionIO.SetData(field,[0,0,0])
            if is_backward:
                Kratos.Expression.LiteralExpressionIO.SetData(field,[0,0,0])
            else:
                Kratos.Expression.NodalPositionExpressionIO.Read(field, Kratos.Configuration.Initial)

            Kratos.Expression.VariableExpressionIO.Write(field, KratosOA.HELMHOLTZ_VECTOR,True)

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.filter_model_part.FullName()}, control variable = {KratosOA.SHAPE.Name()}"
