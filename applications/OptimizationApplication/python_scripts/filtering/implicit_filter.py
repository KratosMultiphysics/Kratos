import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.filtering.filter import Filter
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
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

        if isinstance(filtering_variable, Kratos.DoubleVariable):
            self.filter_variables: 'list[SupportedSensitivityFieldVariableTypes]' = [KratosOA.HELMHOLTZ_SCALAR]
            helmholtz_settings = self.__GetImplicitFilterParameters(filtering_model_part_name, "general_scalar", self.parameters)
        elif isinstance(filtering_variable, Kratos.Array1DVariable3):
            self.filter_variables: 'list[SupportedSensitivityFieldVariableTypes]' = [KratosOA.HELMHOLTZ_VECTOR_X, KratosOA.HELMHOLTZ_VECTOR_Y, KratosOA.HELMHOLTZ_VECTOR_Z]
            if filtering_variable == KratosOA.SHAPE:
                helmholtz_settings = self.__GetImplicitFilterParameters(filtering_model_part_name, "shape", self.parameters)
            else:
                helmholtz_settings = self.__GetImplicitFilterParameters(filtering_model_part_name, "general_vector", self.parameters)
        else:
            raise RuntimeError(f"Unsupported variable = \"{filtering_variable.Name()}\". Only supports DoubleVariable and Array1DVariable3.")

        self.filter_analysis = HelmholtzAnalysis(self.model, helmholtz_settings)

    def Initialize(self) -> None:
        self.filter_analysis.Initialize()
        self.__ApplyBoundaryConditions()

    def Update(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        self.filter_analysis.Finalize()

    def ForwardFilterField(self, control_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        if self.filter_analysis._GetComputingModelPart()[KratosOA.NUMBER_OF_HELMHOLTZ_FILTERS] > 1:
            self.__ReInitializeFilteringModelPart()
        return self.filter_analysis.FilterField(control_field)

    def BackwardFilterField(self, physical_mesh_independent_gradient_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        if self.filter_analysis._GetComputingModelPart()[KratosOA.NUMBER_OF_HELMHOLTZ_FILTERS] > 1:
            self.__ReInitializeFilteringModelPart()
        return self.filter_analysis.FilterField(physical_mesh_independent_gradient_field)

    def BackwardFilterIntegratedField(self, physical_mesh_dependent_gradient_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        if self.filter_analysis._GetComputingModelPart()[KratosOA.NUMBER_OF_HELMHOLTZ_FILTERS] > 1:
            self.__ReInitializeFilteringModelPart()
        return self.filter_analysis.FilterIntegratedField(physical_mesh_dependent_gradient_field)

    def UnfilterField(self, physical_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        if self.filter_analysis._GetComputingModelPart()[KratosOA.NUMBER_OF_HELMHOLTZ_FILTERS] > 1:
            self.__ReInitializeFilteringModelPart()
        return self.UnfilterField(physical_field)

    def __GetImplicitFilterParameters(self, filter_model_part_name: str, filter_type: str, parameters: Kratos.Parameters):
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
        implicit_vector_filter_parameters["solver_settings"]["filter_type"].SetString(filter_type)
        implicit_vector_filter_parameters["solver_settings"]["model_part_name"].SetString(filter_model_part_name)

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
        list(map(lambda filter_var: Kratos.VariableUtils().ApplyFixity(filter_var, False, self.model[self.filtering_model_part_name]), self.filter_variables))

        # now re-apply the boundary conditions
        all_component_boundary_model_parts = KratosOA.OptimizationUtils.GetComponentWiseModelParts(self.model, self.parameters["filtering_boundary_conditions"], len(self.filter_variables))
        list_of_boundary_model_parts: 'list[Kratos.ModelPart]' = []

        for component_index, boundary_model_parts in enumerate(all_component_boundary_model_parts):
            for model_part in boundary_model_parts:
                list_of_boundary_model_parts.append(model_part)
                Kratos.VariableUtils().ApplyFixity(self.filter_variables[component_index], True, model_part.Nodes)

        # make all the boundaries to zero
        for model_part in list_of_boundary_model_parts:
            zero_field =  Kratos.Expression.NodalExpression(model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(zero_field, 0.0)
            for filter_variable in self.filter_variables:
                Kratos.Expression.VariableExpressionIO.Write(zero_field, filter_variable, True)


