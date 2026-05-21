import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as Statistics

from KratosMultiphysics.StatisticsApplication.method_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetMethod

def GetFormattedValue(value, length: int, precision: int) -> str:
    if isinstance(value, int):
        fixed_format = ("{: " + str(length) + ".0f}").format(value)
        if len(fixed_format) <= length:
            return [fixed_format]
        else:
            return [("{: " + str(length) + "." + str(precision) + "e}").format(value)]
    elif isinstance(value, float):
        fixed_format = ("{: " + str(length) + "." + str(precision) + "f}").format(value)
        if len(fixed_format) <= length:
            return [fixed_format]
        else:
            return [("{: " + str(length) + "." + str(precision) + "e}").format(value)]

    elif isinstance(value, Kratos.Array3):
        return [GetFormattedValue(value[0], length, precision)[0], GetFormattedValue(value[1], length, precision)[0], GetFormattedValue(value[2], length, precision)[0]]

def GetItemContainer(item_container_name):
    item_container_types = [
        ["nodal_historical", Statistics.SpatialMethods.Historical],
        ["nodal_non_historical", Statistics.SpatialMethods.NonHistorical.Nodes],
        ["element_non_historical", Statistics.SpatialMethods.NonHistorical.Elements],
        ["condition_non_historical", Statistics.SpatialMethods.NonHistorical.Conditions],
    ]

    item_container_names_list = [item_container_types[i][0] for i in range(len(item_container_types))]
    item_container_types_list = [item_container_types[i][1] for i in range(len(item_container_types))]

    if (item_container_name not in item_container_names_list):
        msg = "Unknown container [ \"container\" = \"" + item_container_name + "\" ]\n"
        msg += "Allowed containers are:\n    "
        msg += "\n    ".join(sorted(item_container_names_list))
        raise Exception(msg)

    return item_container_types_list[item_container_names_list.index(item_container_name)]


def GetVariableHeaders(norm_type, variable_name):
    if (norm_type == "none"):
        variable_type = Kratos.KratosGlobals.GetVariableType(variable_name)
        if (variable_type == "Double"):
            return [variable_name]
        elif (variable_type == "Array" and Kratos.KratosGlobals.GetVariableType(variable_name + "_X") == "Double"):
            return [variable_name + "_X", variable_name + "_Y", variable_name + "_Z"]
        else:
            raise Exception("Unsupported variable type " + variable_type)
    else:
        return [variable_name]

class SpatialMethodOutput:
    def __init__(self, method_name, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        self.method_name = method_name
        self.container_type = container_type
        self.norm_type = norm_type
        self.variable_name = variable_name
        self.method_settings = method_settings
        self.variable = Kratos.KratosGlobals.GetVariable(self.variable_name)

        self.container = GetItemContainer(self.container_type)
        self.norm_type_container = GetNormTypeContainer(self.container, norm_type)
        method = GetMethod(self.norm_type_container, self.method_name)

        if method_settings.IsEquivalentTo(Kratos.Parameters("""{}""")):
            if norm_type == "none":
                self.norm_type_container = self.container.ValueMethods
                self.method = lambda model_part: method(model_part, self.variable)
            else:
                self.norm_type_container = self.container.NormMethods
                self.method = lambda model_part: method(model_part, self.variable, norm_type)
        else:
            if norm_type == "none":
                self.norm_type_container = self.container.ValueMethods
                self.method = lambda model_part: method(model_part, self.variable, method_settings)
            else:
                self.norm_type_container = self.container.NormMethods
                self.method = lambda model_part: method(model_part, self.variable, norm_type, method_settings)

    def Evaluate(self, model_part):
        self.data = self.method(model_part)

    def GetValues(self, length, precision):
        if isinstance(self.data, float):
            return GetFormattedValue(self.data, length, precision)
        else:
            values = []
            for v in self.data:
                values.extend(GetFormattedValue(v, length, precision))
            return values

    def GetValueLengths(self, value_length):
        return [value_length]

    def __str__(self):
        return self.__class__.__name__ + "(container: {:s}, norm_type: {:s}, variable: {:s})".format(self.container_type, self.norm_type, self.variable_name)


class SpatialSumOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        super().__init__("sum", container_type, norm_type, variable_name)

    def GetHeaders(self):
        return ["sum"]


class SpatialMeanOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        super().__init__("mean", container_type, norm_type, variable_name)

    def GetHeaders(self):
        return ["mean"]


class SpatialRootMeanSquareOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        super().__init__("rootmeansquare", container_type, norm_type, variable_name)

    def GetHeaders(self):
        return ["rootmeansquare"]


class SpatialMinOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        super().__init__("min", container_type, norm_type, variable_name)

    def GetHeaders(self):
        return ["min_value", "min_id"]

    def GetValueLengths(self, value_length):
        return [value_length, value_length]


class SpatialMaxOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        super().__init__("max", container_type, norm_type, variable_name)

    def GetHeaders(self):
        return ["max_value", "max_id"]

    def GetValueLengths(self, value_length):
        return [value_length, value_length]


class SpatialMedianOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        super().__init__("median", container_type, norm_type, variable_name)

    def GetHeaders(self):
        return ["median"]


class SpatialVarianceOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings = Kratos.Parameters("""{}""")):
        super().__init__("variance", container_type, norm_type, variable_name)

    def GetHeaders(self):
        return ["mean", "variance"]

    def GetValueLengths(self, value_length):
        return [value_length, value_length]


class SpatialDistributionOutput(SpatialMethodOutput):
    def __init__(self, container_type, norm_type, variable_name, method_settings):
        super().__init__("distribution", container_type, norm_type, variable_name, method_settings)
        default_parameters = Kratos.Parameters("""{ "number_of_value_groups": 10 }""")
        self.method_settings.AddMissingParameters(default_parameters)

        self.write_min_value = True
        if self.method_settings.Has("min_value") and self.method_settings["min_value"].IsDouble():
            self.write_min_value = False

        self.write_max_value = True
        if self.method_settings.Has("max_value") and self.method_settings["max_value"].IsDouble():
            self.write_max_value = False

    def Evaluate(self, *args):
        super().Evaluate(*args)
        (self.min_value, self.max_value, self.group_upper_values, self.group_histogram, self.group_percentage_distribution, self.group_means, self.group_variances) = self.data

    def GetHeaders(self):
        number_of_groups = self.method_settings["number_of_value_groups"].GetInt()
        headers = []
        if self.write_min_value:
            headers.append("min")
        if self.write_max_value:
            headers.append("max")

        headers.append("group_below_min {mean|variance}")
        headers.extend(["group_{:d}".format(i+1) + " {mean|variance}" for i in range(number_of_groups)])
        headers.append("group_above_max {mean|variance}")
        return headers

    def GetValues(self, length, precision):
        values = []
        if self.write_min_value:
            values.append(GetFormattedValue(self.min_value, length, precision)[0])
        if self.write_max_value:
            values.append(GetFormattedValue(self.max_value, length, precision)[0])

        def get_formatted_group_value(hist_v, percentage_v, upper_v, mean_v, variance_v):
            s = "{:s} [#] = {:s} [%] < {:s}".format(GetFormattedValue(hist_v, length, precision)[0], GetFormattedValue(percentage_v * 100.0, length, precision)[0], GetFormattedValue(upper_v, length, precision)[0])
            s += " { " + GetFormattedValue(mean_v, length, precision)[0] + " | " + GetFormattedValue(variance_v, length, precision)[0] + " }"
            return s
        values.extend([get_formatted_group_value(hist_v, percentage_v, upper_v, mean_v, variance_v) for hist_v, percentage_v, upper_v, mean_v, variance_v in zip(self.group_histogram, self.group_percentage_distribution, self.group_upper_values, self.group_means, self.group_variances)])
        values[-1] = values[-1].replace(" <", ">=")
        return values

    def GetValueLengths(self, value_length):
        values_per_header = []
        if self.write_min_value:
            values_per_header.append(value_length)
        if self.write_max_value:
            values_per_header.append(value_length)

        number_of_groups = self.method_settings["number_of_value_groups"].GetInt()
        for _ in range(number_of_groups + 2):
            values_per_header.append(22 + value_length * 5)
        return values_per_header


