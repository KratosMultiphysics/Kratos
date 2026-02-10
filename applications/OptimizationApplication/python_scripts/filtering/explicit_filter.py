import operator
from itertools import accumulate

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.filtering.filter import Filter
from KratosMultiphysics.OptimizationApplication.filtering.filter_utils import FilterRadiusFactory
from KratosMultiphysics.OptimizationApplication.filtering.filter_utils import ExplicitFilterDampingFactory
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

def Factory(model: Kratos.Model, filtering_model_part_name: str, filtering_variable: SupportedSensitivityFieldVariableTypes, data_location: Kratos.Globals.DataLocation, parameters: Kratos.Parameters) -> Filter:
    return ExplicitFilter(model, filtering_model_part_name, filtering_variable, data_location, parameters)

class ExplicitFilter(Filter):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "filter_type"               : "explicit_filter",
            "filter_function_type"      : "linear",
            "max_items_in_bucket"       : 10,
            "echo_level"                : 0,
            "filter_radius_settings":{
                "filter_radius_type": "constant",
                "filter_radius"     : 0.2
            },
            "filtering_boundary_conditions": {
                "damping_type"              : "nearest_entity",
                "damping_function_type"     : "cosine",
                "damped_model_part_settings": {}
            }
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
            self.filter_variable_shape = []
        elif isinstance(filtering_variable, Kratos.Array1DVariable3):
            self.filter_variable_shape = [3]
        else:
            raise RuntimeError(f"Unsupported variable = \"{filtering_variable.Name()}\". Only supports DoubleVariable and Array1DVariable3.")

    def Initialize(self) -> None:
        # get the model part
        self.model_part = self.model.GetModelPart(self.filtering_model_part_name)

        # we create the filter in the initialize to support filter to work
        # on sub model parts.
        filter_function_type = self.parameters["filter_function_type"].GetString()
        max_items_in_bucket = self.parameters["max_items_in_bucket"].GetInt()
        echo_level = self.parameters["echo_level"].GetInt()
        if self.data_location in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
            self.filter_utils = KratosOA.NodeExplicitFilterUtils(self.model_part, filter_function_type, max_items_in_bucket, echo_level)
        elif self.data_location == Kratos.Globals.DataLocation.Condition:
            self.filter_utils = KratosOA.ConditionExplicitFilterUtils(self.model_part, filter_function_type, max_items_in_bucket, echo_level)
        elif self.data_location == Kratos.Globals.DataLocation.Element:
            self.filter_utils = KratosOA.ElementExplicitFilterUtils(self.model_part, filter_function_type, max_items_in_bucket, echo_level)

        self.Update()

    def Update(self) -> None:
        # now set the filter radius. Can be changed in future to support adaptive radius methods.
        filter_radius = FilterRadiusFactory(self.model_part, self.data_location, self.parameters["filter_radius_settings"])
        self.filter_utils.SetRadius(filter_radius)
        self.GetComponentDataView().GetUnBufferedData().SetValue("filter_radius", filter_radius.Clone(), overwrite=True)

        # now set the damping
        stride = max(enumerate(accumulate(self.filter_variable_shape, operator.mul, initial=1)))[1]
        self.damping = ExplicitFilterDampingFactory(self.model, self.data_location, stride, self.parameters["filtering_boundary_conditions"])
        self.damping.SetRadius(filter_radius)
        self.filter_utils.SetDamping(self.damping)

        # initialize the filter
        self.damping.Update()
        self.filter_utils.Update()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def ForwardFilterField(self, control_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return self.filter_utils.ForwardFilterField(control_field)

    def BackwardFilterField(self, physical_mesh_independent_gradient_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return self.filter_utils.BackwardFilterField(physical_mesh_independent_gradient_field)

    def BackwardFilterIntegratedField(self, physical_mesh_dependent_gradient_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return self.filter_utils.BackwardFilterIntegratedField(physical_mesh_dependent_gradient_field)

    def UnfilterField(self, filtered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
         # WARNING: In general explicit VM doesn't need unfiltered field because it works with changes. Hence, it returns zeros and optimization runs correctly.
        return Kratos.Expression.Utils.Collapse(filtered_field * 0.0)

    def GetBoundaryConditions(self) -> 'list[list[Kratos.ModelPart]]':
        return self.damping.GetDampedModelParts()