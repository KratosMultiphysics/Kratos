import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.filtering.filter import Filter
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis

def Factory(model: Kratos.Model, filtering_model_part_name: str, filtering_variable: SupportedSensitivityFieldVariableTypes, data_location: Kratos.Globals.DataLocation, parameters: Kratos.Parameters) -> Filter:
    return ImplicitFilter(model, filtering_model_part_name, filtering_variable, data_location, parameters)

class ImplicitFilter(Filter):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "filter_type"                  : "implicit_filter",
            "filter_radius"                : 0.2,
            "filter_variable_type"         : "auto",
            "linear_solver_settings"       : {},
            "echo_level"                   : 0,
            "filtering_boundary_conditions": {}
        }""")

    def __init__(self, model: Kratos.Model, filtering_model_part_name: str, filtering_variable: SupportedSensitivityFieldVariableTypes, data_location: Kratos.Globals.DataLocation, parameters: Kratos.Parameters) -> None:
        super().__init__()

        self.model = model
        self.filtering_model_part_name = filtering_model_part_name
        self.filtering_variable = filtering_variable
        self.data_location = data_location

        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        supported_data_locations = [
            Kratos.Globals.DataLocation.NodeHistorical,
            Kratos.Globals.DataLocation.NodeNonHistorical,
            Kratos.Globals.DataLocation.Condition,
            Kratos.Globals.DataLocation.Element
        ]
        if data_location not in supported_data_locations:
            raise RuntimeError(f"Unsupported data location = \"{data_location.name}\" requested. Followings are supported:\n\t" + "\n\t".join([v.name for v in supported_data_locations]))

        allowed_filter_variable_types = ["auto", "general_scalar", "general_vector", "shape"]
        filter_variable_type = self.parameters["filter_variable_type"].GetString()
        if not filter_variable_type in allowed_filter_variable_types:
            raise RuntimeError(f"Unsupported filter_variable_type = \"{filter_variable_type}\". Followings are supported:\n\t" + "\n\t".join(allowed_filter_variable_types))

        if filter_variable_type == "auto":
            if isinstance(filtering_variable, Kratos.DoubleVariable):
                self.filter_type = "general_scalar"
            elif isinstance(filtering_variable, Kratos.Array1DVariable3):
                if filtering_variable == KratosOA.SHAPE:
                    self.filter_type = "shape"
                else:
                    self.filter_type = "general_vector"
            else:
                raise RuntimeError(f"Unsupported variable type [ variable name = {filtering_variable.Name()}]")
        else:
            self.filter_type = filter_variable_type

        if self.filter_type == "general_scalar" and not isinstance(filtering_variable, Kratos.DoubleVariable):
            raise RuntimeError(f"\"general_scalar\" only supports scalar variables [ variable name = {filtering_variable.Name()} ].")
        if self.filter_type == "general_vector" and not isinstance(filtering_variable, Kratos.Array1DVariable3):
            raise RuntimeError(f"\"general_vector\" only supports array3 variables [ variable name = {filtering_variable.Name()} ].")
        if self.filter_type == "shape" and filtering_variable != KratosOA.SHAPE:
            raise RuntimeError(f"\"shape\" only supports SHAPE variable [ variable name = {filtering_variable.Name()} ].")

        helmholtz_settings = self.__GetImplicitFilterParameters(filtering_model_part_name, self.parameters)

        self.filter_analysis = HelmholtzAnalysis(self.model, helmholtz_settings)

        solving_var = self.filter_analysis._GetSolver().GetSolvingVariable()
        if isinstance(solving_var, Kratos.DoubleVariable):
            self.filter_variables = [solving_var]
        elif isinstance(solving_var, Kratos.Array1DVariable3):
            self.filter_variables = list(map(lambda x: Kratos.KratosGlobals.GetVariable(x), [f"{solving_var.Name()}_{suffix}" for suffix in ["X", "Y", "Z"]]))
        else:
            raise RuntimeError(f"Unsupported solving variable [ solving variable name = {solving_var.Name()} ].")

        self.damped_model_parts: 'typing.Optional[list[list[Kratos.ModelPart]]]' = None

    def Initialize(self) -> None:
        self.filter_analysis.Initialize()
        self.__ApplyBoundaryConditions()

    def Update(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        self.filter_analysis.Finalize()

    def ForwardFilterField(self, control_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        if self.filter_analysis._GetSolver().GetOriginModelPart()[KratosOA.NUMBER_OF_SOLVERS_USING_NODES] > 1:
            self.__ReInitializeFilteringModelPart()
        return self.filter_analysis.FilterField(control_field)

    def BackwardFilterField(self, physical_mesh_independent_gradient_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        if self.filter_analysis._GetSolver().GetOriginModelPart()[KratosOA.NUMBER_OF_SOLVERS_USING_NODES] > 1:
            self.__ReInitializeFilteringModelPart()
        return self.filter_analysis.FilterField(physical_mesh_independent_gradient_field)

    def BackwardFilterIntegratedField(self, physical_mesh_dependent_gradient_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        if self.filter_analysis._GetSolver().GetOriginModelPart()[KratosOA.NUMBER_OF_SOLVERS_USING_NODES] > 1:
            self.__ReInitializeFilteringModelPart()
        return self.filter_analysis.FilterIntegratedField(physical_mesh_dependent_gradient_field)

    def UnfilterField(self, physical_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        if self.filter_analysis._GetSolver().GetOriginModelPart()[KratosOA.NUMBER_OF_SOLVERS_USING_NODES] > 1:
            self.__ReInitializeFilteringModelPart()
        # now we get the current values
        current_values_ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.filter_analysis._GetComputingModelPart().Nodes, self.filter_analysis._GetSolver().GetSolvingVariable())
        current_values_ta.CollectData()

        # now apply the physical field
        self.filter_analysis.AssignTensorDataToNodalSolution(physical_field)
        unfiltered_field = self.filter_analysis.UnFilterField(physical_field)

        # now apply back the original BCs
        current_values_ta.StoreData()
        return unfiltered_field

    def GetBoundaryConditions(self) -> 'list[list[Kratos.ModelPart]]':
        return self.damped_model_parts

    def __GetImplicitFilterParameters(self, filter_variable_typel_part_name: str, parameters: Kratos.Parameters):
        implicit_vector_filter_parameters = Kratos.Parameters("""
        {
            "solver_settings" : {
                "domain_size"          : 3,
                "echo_level"           : 0,
                "filter_type"          : "general_vector",
                "model_part_name"      : "",
                "filter_radius"        : 0.0,
                "model_import_settings": {
                    "input_type"     : "use_input_model_part"
                }
            },
            "problem_data": {
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 1.0,
                "parallel_type" : "OpenMP"
            }
        }""")

        implicit_vector_filter_parameters["solver_settings"]["echo_level"].SetInt(parameters["echo_level"].GetInt())
        implicit_vector_filter_parameters["solver_settings"]["filter_radius"].SetDouble(parameters["filter_radius"].GetDouble())
        implicit_vector_filter_parameters["solver_settings"]["filter_type"].SetString(self.filter_type)
        implicit_vector_filter_parameters["solver_settings"]["model_part_name"].SetString(filter_variable_typel_part_name)

        linear_solver_settings = parameters["linear_solver_settings"]
        if linear_solver_settings.Has("solver_type"):
            implicit_vector_filter_parameters["solver_settings"].AddValue("linear_solver_settings", linear_solver_settings)

        return implicit_vector_filter_parameters

    def __ReInitializeFilteringModelPart(self) -> None:
        self.filter_analysis.InitializeFilterModelPart()
        self.__ApplyBoundaryConditions()

    def __ApplyBoundaryConditions(self) -> None:
        # first free the domain, because if the same type of filter used, that means same nodes
        # were used before with different boundary conditions. Hence it is important to free all dofs
        # and apply the boundary conditions again when ever filter is used.
        # This is because, every implicit filter shares the same nodes
        list(map(lambda filter_var: Kratos.VariableUtils().ApplyFixity(filter_var, False, self.filter_analysis._GetSolver().GetComputingModelPart().Nodes), self.filter_variables))

        # now re-apply the boundary conditions
        self.damped_model_parts = KratosOA.OptimizationUtils.GetComponentWiseModelParts(self.model, self.parameters["filtering_boundary_conditions"])
        if len(self.filter_variables) != len(self.damped_model_parts) and len(self.damped_model_parts) != 0:
            raise RuntimeError("Number of components mismatch.")
        list_of_boundary_model_parts: 'list[Kratos.ModelPart]' = []

        for component_index, boundary_model_parts in enumerate(self.damped_model_parts):
            for model_part in boundary_model_parts:
                list_of_boundary_model_parts.append(model_part)
                Kratos.VariableUtils().ApplyFixity(self.filter_variables[component_index], True, model_part.Nodes)

        # make all the boundaries to zero
        for model_part in list_of_boundary_model_parts:
            for filter_variable in self.filter_variables:
                zero_field = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(model_part.Nodes, filter_variable)
                zero_field.data[:] = 0.0
                zero_field.StoreData()

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Applied filter boundary conditions used in \"{self.GetComponentDataView().GetComponentName()}\".")


