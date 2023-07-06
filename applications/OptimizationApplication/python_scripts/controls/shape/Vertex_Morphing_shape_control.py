import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"VertexMorphingShapeControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"VertexMorphingShapeControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return VertexMorphingShapeControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

def GetImplicitFilterParameters(model_part: Kratos.ModelPart, parameters: Kratos.Parameters):
    model_part_name = model_part.FullName()
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
        implicit_vector_filter_parameters["solver_settings"]["filter_type"] = "bulk_surface_shape"

    implicit_vector_filter_parameters["solver_settings"].AddDouble("filter_radius",filter_radius)
    implicit_vector_filter_parameters["solver_settings"].AddString("model_part_name",model_part_name)
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
                "type" : "surface_implicit",
                "radius": 0.000000000001,
                "linear_solver_settings" : {}
            },
            "output_all_fields": false,
            "fixed_model_parts": {},
            "utilities": []
        }""")

        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.parameters["filter_settings"].ValidateAndAssignDefaults(default_settings["filter_settings"])

        self.output_all_fields = self.parameters["output_all_fields"].GetBool()
        self.filter_type = self.parameters["filter_settings"]["type"].GetString()

        self.supported_filter_types = ["surface_implicit","bulk_surface_implicit","explicit"]
        if not self.filter_type in self.supported_filter_types:
            raise RuntimeError(f"The specified filter type \"{self.filter_type}\" is not supported. Followings are the variables:\n\t" + "\n\t".join([k for k in self.supported_filter_types]))

        controlled_model_parts = [model[model_part_name] for model_part_name in parameters["controlled_model_part_names"].GetStringArray()]
        if len(controlled_model_parts) == 0:
            raise RuntimeError(f"No model parts were provided for VertexMorphingShapeControl. [ control name = \"{self.GetName()}\"]")

        self.model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.UNION, f"control_{self.GetName()}", controlled_model_parts, False)

        self.fixed_model_parts: 'dict[str, list]' = {}
        for model_part_name, direction_params in self.parameters["fixed_model_parts"].items():
            bool_list = []
            for i in range(3):
                bool_list.append(direction_params[i].GetBool())
            self.fixed_model_parts[model_part_name] = bool_list

        self.filter = None
        if "implicit" in self.filter_type:
            helmholtz_settings = GetImplicitFilterParameters(self.model_part, self.parameters)
            self.filter = HelmholtzAnalysis(self.model, helmholtz_settings)

        self.is_initialized = False

    def Initialize(self) -> None:
        ModelPartUtilities.ExecuteOperationOnModelPart(self.model_part)

        # initiliaze the filter
        self.filter.Initialize()
        physical_shape_field = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(physical_shape_field, Kratos.Configuration.Initial)
        self._SetFixedModelPartValues(False)
        self.control_field = self.filter.UnFilterField(physical_shape_field)

        # now update
        self._UpdateAndOutputFields()
        self.is_initialized = True

    def Check(self) -> None:
        return self.filter.Check()

    def Finalize(self) -> None:
        return self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, [0,0,0])
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field

    def GePhysicalField(self) -> ContainerExpressionTypes:
        physical_shape_field = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(physical_shape_field, Kratos.Configuration.Initial)
        return physical_shape_field

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        with TimeLogger("VertexMorphingShapeControl::MapGradient", None, "Finished",False):
            keys = physical_gradient_variable_container_expression_map.keys()
            if len(keys) != 1:
                raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
            if KratosOA.SHAPE not in keys:
                raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {KratosOA.SHAPE.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

            physical_gradient = physical_gradient_variable_container_expression_map[KratosOA.SHAPE]
            if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
                raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

            # now filter the field
            filtered_gradient = self.filter.FilterIntegratedField(physical_gradient)

            return filtered_gradient

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")

        if KratosOA.ExpressionUtils.NormL2(self.control_field - new_control_field) > 1e-15:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}.",False):
                # update the control SHAPE field
                self.control_field = new_control_field
                # now update the physical field
                self._UpdateAndOutputFields()
                return True
        return False

    def _UpdateAndOutputFields(self) -> None:

        # apply boundary conditions for the forward filtering
        self._SetFixedModelPartValues(False)
        # filter the control field
        new_shape = self.filter.FilterField(self.control_field)
        # apply boundary conditions for the filtering of the sensitivites
        self._SetFixedModelPartValues(True)

        # compute shape update
        shape_update = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(shape_update, Kratos.Configuration.Initial)
        shape_update -= new_shape

        # now update shape
        Kratos.Expression.NodalPositionExpressionIO.Write(new_shape, Kratos.Configuration.Initial)
        Kratos.Expression.NodalPositionExpressionIO.Write(new_shape, Kratos.Configuration.Current)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("shape_update", shape_update.Clone(),overwrite=self.is_initialized)
        if self.output_all_fields:
            un_buffered_data.SetValue("control_shape", self.control_field.Clone(),overwrite=self.is_initialized)

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
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {KratosOA.SHAPE.Name()}"
