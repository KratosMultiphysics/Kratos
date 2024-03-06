import typing
import operator
from itertools import accumulate

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.filtering.filter import Filter
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

def Factory(model: Kratos.Model, filtering_model_part_name: str, filtering_variable: SupportedSensitivityFieldVariableTypes, data_location: Kratos.Globals.DataLocation, parameters: Kratos.Parameters) -> Filter:
    return ExplicitVertexMorphingFilter(model, filtering_model_part_name, filtering_variable, data_location, parameters)

class ExplicitVertexMorphingFilter(Filter):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "filter_type"               : "explicit_vertex_morphing",
            "filter_function_type"      : "linear",
            "max_nodes_in_filter_radius": 100000,
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
        max_nodes_in_filter_radius = self.parameters["max_nodes_in_filter_radius"].GetInt()
        if self.data_location in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
            self.filter = KratosOA.NodalExplicitFilter(self.model_part, filter_function_type, max_nodes_in_filter_radius)
        elif self.data_location == Kratos.Globals.DataLocation.Condition:
            self.filter = KratosOA.ConditionExplicitFilter(self.model_part, filter_function_type, max_nodes_in_filter_radius)
        elif self.data_location == Kratos.Globals.DataLocation.Element:
            self.filter = KratosOA.ElementExplicitFilter(self.model_part, filter_function_type, max_nodes_in_filter_radius)

        self.Update()

    def Update(self) -> None:
        # now set the filter radius. Can be changed in future to support adaptive radius methods.
        filter_radius = self._GetFilterRadiusExpression(self.parameters["filter_radius_settings"])
        self.filter.SetFilterRadius(filter_radius)
        self.GetComponentDataView().GetUnBufferedData().SetValue("filter_radius", filter_radius.Clone(), overwrite=True)

        # now set the damping
        damping_coefficients = self._GetDampingCoefficients(filter_radius, self.parameters["filtering_boundary_conditions"])
        self.filter.SetDampingCoefficients(damping_coefficients)
        self.GetComponentDataView().GetUnBufferedData().SetValue("damping_coefficients", damping_coefficients.Clone(), overwrite=True)

        # initialize the filter
        self.filter.Update()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def FilterField(self, unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return self.filter.FilterField(unfiltered_field)

    def FilterIntegratedField(self, unfiltered_integrated_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return self.filter.FilterIntegratedField(unfiltered_integrated_field)

    def _GetFilterRadiusExpression(self, filter_radius_settings: Kratos.Parameters) -> ContainerExpressionTypes:
        if not filter_radius_settings.Has("filter_radius_type"):
            raise RuntimeError(f"\"filter_radius_type\" not specified in the following settings:\n{filter_radius_settings}")

        filter_radius_type = filter_radius_settings["filter_radius_type"].GetString()
        if filter_radius_type == "constant":
            defaults = Kratos.Parameters("""{
                "filter_radius_type": "constant",
                "filter_radius"     : 0.2
            }""")
            filter_radius_settings.ValidateAndAssignDefaults(defaults)
            filter_radius = self._GetContainerExpressionType()(self.model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, filter_radius_settings["filter_radius"].GetDouble())
            return filter_radius
        else:
            raise RuntimeError(f"Unsupported filter_radius_type = \"{filter_radius_type}\".")

    def _GetDampingCoefficients(self, damping_radius:ContainerExpressionTypes, filter_boundary_condition_settings: Kratos.Parameters) -> ContainerExpressionTypes:
        if not filter_boundary_condition_settings.Has("damping_type"):
            raise RuntimeError(f"\"damping_type\" not specified in the following settings:\n{filter_boundary_condition_settings}")

        damping_type = filter_boundary_condition_settings["damping_type"].GetString()
        if damping_type == "nearest_entity":
            defaults = Kratos.Parameters("""{
                "damping_type"              : "nearest_entity",
                "damping_function_type"     : "cosine",
                "bucket_size"               : 100,
                "damped_model_part_settings": {}
            }""")
            filter_boundary_condition_settings.ValidateAndAssignDefaults(defaults)

            number_of_components = max(enumerate(accumulate(self.filter_variable_shape, operator.mul, initial=1)))[1]
            damped_model_parts = KratosOA.FilterUtils.GetComponentWiseDampedModelParts(self.model, filter_boundary_condition_settings["damped_model_part_settings"], number_of_components)
            return KratosOA.FilterUtils.ComputeDampingCoefficientsBasedOnNearestEntity(damping_radius, damped_model_parts, self.filter_variable_shape, filter_boundary_condition_settings["damping_function_type"].GetString(), filter_boundary_condition_settings["bucket_size"].GetInt())
        else:
            raise RuntimeError(f"Unsupported damping_type = \"{damping_type}\".")

    def _GetContainerExpressionType(self) -> typing.Type[ContainerExpressionTypes]:
        if self.data_location in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
            return Kratos.Expression.NodalExpression
        elif self.data_location == Kratos.Globals.DataLocation.Condition:
            return Kratos.Expression.ConditionExpression
        elif self.data_location == Kratos.Globals.DataLocation.Element:
            return Kratos.Expression.ElementExpression



