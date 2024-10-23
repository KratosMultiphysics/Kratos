from abc import ABC, abstractmethod
from typing import Any, Union
from math import log10

import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStat

def GetListOfValues(value: Any) -> 'Union[list[int],list[float]]':
    if isinstance(value, list):
        return value
    elif isinstance(value, int) or isinstance(value, float):
        return [value]
    elif isinstance(value, Kratos.Array3) or isinstance(value, Kratos.Array4) or isinstance(value, Kratos.Array6) or isinstance(value, Kratos.Array9) or isinstance(value, Kratos.Vector):
        return [v for v in value]
    elif isinstance(value, Kratos.Matrix):
        values: 'list[float]' = []
        for i in range(value.Size1()):
            for j in range(value.Size2()):
                values.append(value[i, j])
        return values
    else:
        raise RuntimeError(f"Unsupported value type provided [ value = {value} ]. Supported value types are:\n\tint\n\tfloat\n\tArray3\n\tArray4\n\tArray6\n\tArray9\n\tVector\n\tMatrix")

def GetListOfValueHeaderSuffixes(value: Any) -> 'list[str]':
    if isinstance(value, int) or isinstance(value, float):
        return [""]
    elif isinstance(value, Kratos.Array3):
        return ["_X", "_Y", "_Z"]
    elif isinstance(value, Kratos.Array4) or isinstance(value, Kratos.Array6) or isinstance(value, Kratos.Array9) or isinstance(value, Kratos.Vector):
        return [f"_[{i+1}]" for i, _ in enumerate(value)]
    elif isinstance(value, Kratos.Matrix):
        values: 'list[str]' = []
        for i in range(value.Size1()):
            for j in range(value.Size2()):
                values.append(f"_[{i+1,j+1}]")
        return values
    else:
        raise RuntimeError(f"Unsupported value type provided [ value = {value} ]. Supported value types are:\n\tint\n\tfloat\n\tArray3\n\tArray4\n\tArray6\n\tArray9\n\tVector\n\tMatrix")

def GetDataRetrievalMethods(container_location: Kratos.Globals.DataLocation) -> Any:
    if container_location == Kratos.Globals.DataLocation.NodeHistorical:
        return lambda x: x.Nodes, lambda x, y: x.GetSolutionStepValue(y)
    elif container_location == Kratos.Globals.DataLocation.NodeNonHistorical:
        return lambda x: x.Nodes, lambda x, y: x.GetValue(y)
    elif container_location == Kratos.Globals.DataLocation.Condition:
        return lambda x: x.Conditions, lambda x, y: x.GetValue(y)
    elif container_location == Kratos.Globals.DataLocation.Element:
        return lambda x: x.Elements, lambda x, y: x.GetValue(y)
    else:
        raise RuntimeError("Unsupported container location.")

def GetValueForVariable(variable: Any) -> Any:
    if isinstance(variable, Kratos.IntegerVariable):
        return 0
    elif isinstance(variable, Kratos.DoubleVariable):
        return 0.0
    elif isinstance(variable, Kratos.Array1DVariable3):
        return Kratos.Array3(0.0)
    elif isinstance(variable, Kratos.Array1DVariable4):
        return Kratos.Array4(0.0)
    elif isinstance(variable, Kratos.Array1DVariable6):
        return Kratos.Array6(0.0)
    elif isinstance(variable, Kratos.Array1DVariable9):
        return Kratos.Array9(0.0)
    elif isinstance(variable, Kratos.Vector):
        return Kratos.Vector()
    elif isinstance(variable, Kratos.Matrix):
        return Kratos.Matrix()
    else:
        raise RuntimeError(f"Unsupported variable [ variable = {variable} ].")

def GetValueExponentLength(value) -> int:
    abs_value = abs(value)
    if abs_value >= 1e+100:
        return int(log10(log10(abs_value))) + 1
    else:
        return 2

def GetFormattedInt(value: int, precision: int) -> str:
    fixed_scientific_str = ("{: ." + str(precision - GetValueExponentLength(value)) + "e}").format(value)
    fixed_str = ("{:> " + str(len(fixed_scientific_str)) + "d}").format(value)
    if len(fixed_str) <= len(fixed_scientific_str):
        return fixed_str
    else:
        return fixed_scientific_str

def GetFormattedFloat(value: float, precision: int) -> str:
    fixed_scientific_str = ("{: ." + str(precision - GetValueExponentLength(value)) + "e}").format(value)
    fixed_str = ("{: ." + str(precision) + "f}").format(value)
    if len(fixed_str) <= len(fixed_scientific_str):
        return fixed_str
    else:
        return fixed_scientific_str

def GetLengthAdjustedValue(value: 'Union[int, float]', length: int, precision: int) -> str:
    if isinstance(value, int):
        v_str = GetFormattedInt(value, precision)
    elif isinstance(value, float):
        v_str = GetFormattedFloat(value, precision)
    else:
        raise RuntimeError("Unsupported data type")

    if len(v_str) < length:
        v_str = ("{:>" + str(length) + "s}").format(v_str)
    return v_str

def GetTaggedHeader(header_details: 'dict[str, bool]', location_name: str, variable_name: str, norm: Any, suffix: str, method_name: str) -> str:
    header = ""
    if header_details["location"]:
        header += f"@{location_name} "
    if header_details["norm"] and norm is not None:
        header += f"{norm} "
    if header_details["variable"]:
        header += f"{variable_name}{suffix} "
    else:
        if suffix != "":
            header += f"{suffix[1:]} "
    if header_details["method"]:
        header += f"{method_name} <TAG> "
    else:
        header += "<TAG> "
    header = header.strip()
    if header.find(",") != -1:
        header = f"\"{header}\""
    return header

def GetShapeSynchronizedValue(model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation):
    container_retrieval_method, entity_data_retireval_method = GetDataRetrievalMethods(container_location)
    dummy_value = GetValueForVariable(variable)
    for entity in container_retrieval_method(model_part):
        dummy_value = entity_data_retireval_method(entity, variable)
        break
    model_part.GetCommunicator().GetDataCommunicator().SynchronizeShape(dummy_value)
    return dummy_value

class SpatialStatisticsOperation(ABC):
    def __init__(self, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any) -> None:
        self.model_part = model_part
        self.variable = variable
        self.container_location = container_location
        self.norm = norm

    def Check(self) -> None:
        if self.container_location == Kratos.Globals.DataLocation.NodeHistorical:
            if not self.model_part.HasNodalSolutionStepVariable(self.variable):
                raise RuntimeError(f"The statistics output operation with {self.variable.Name()} requires the {self.model_part.FullName()} to have the variable in the nodal solution step variables list.")

    def GetNormInfo(self) -> str:
        return str(self.norm)

    def GetVarianbleInfo(self) -> str:
        return self.variable.Name()

    def GetContainerInfo(self) -> str:
        return self.container_location.name

    def _GetSuffixedTaggedHeaders(self, header_details: 'dict[str, bool]', method_name: str) -> 'list[str]':
        if self.norm is None:
            dummy_value = GetShapeSynchronizedValue(self.model_part, self.variable, self.container_location)
        else:
            dummy_value = 0.0

        list_of_tagged_headers = []
        for suffix in GetListOfValueHeaderSuffixes(dummy_value):
            header = GetTaggedHeader(header_details, self.container_location.name, self.variable.Name(), self.norm, suffix, method_name)
            list_of_tagged_headers.append(header)
        return list_of_tagged_headers

    def __str__(self) -> str:
        return f"{self.GetVarianbleInfo()} - {self.GetNormInfo()} - {self.GetContainerInfo()} - {self.GetMethodInfo()}"

    @abstractmethod
    def GetMethodInfo(self) -> str:
        pass

    @abstractmethod
    def GetValueString(self) -> str:
        pass

    @abstractmethod
    def GetHeadersString(self, header_details: 'dict[str, bool]') -> 'str':
        pass

class SpatialStatisticsValueOperation(SpatialStatisticsOperation):
    def __init__(self, method: Any, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, precision: int) -> None:
        super().__init__(model_part, variable, container_location, norm)
        self.method = method
        self.precision = precision
        self.header_lengths: 'list[int]' = []

    def GetMethodInfo(self) -> str:
        return self.method.__name__.lower()

    def GetValueString(self) -> str:
        if self.norm is not None:
            data = self.method(self.model_part, self.variable, self.container_location, self.norm)
        else:
            data = self.method(self.model_part, self.variable, self.container_location)

        list_of_values = GetListOfValues(data)

        list_of_str_values: 'list[str]' = []
        for i, v in enumerate(list_of_values):
            v_str = GetLengthAdjustedValue(v, self.header_lengths[i], self.precision)
            list_of_str_values.append(v_str)
        return ", ".join(list_of_str_values)

    def GetHeadersString(self, header_details: 'dict[str, bool]') -> 'str':
        value_length = GetFormattedFloat(1e+100, self.precision)

        list_of_headers = []
        for header in self._GetSuffixedTaggedHeaders(header_details, self.GetMethodInfo()):
            header = header.replace("<TAG>", "").strip()
            if len(header) < len(value_length):
                header = ("{:>" + str(len(value_length)) + "s}").format(header)
            list_of_headers.append(header)
            self.header_lengths.append(len(header))

        return ", ".join(list_of_headers)

class SpatialStatisticsValueIndexPairOperation(SpatialStatisticsOperation):
    def __init__(self, method: Any, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, precision: int) -> None:
        super().__init__(model_part, variable, container_location, norm)
        self.method = method
        self.precision = precision
        self.header_lengths: 'list[int]' = []

    def GetMethodInfo(self) -> str:
        return self.method.__name__.lower()

    def GetValueString(self) -> str:
        if self.norm is not None:
            values, indices = self.method(self.model_part, self.variable, self.container_location, self.norm)
        else:
            values, indices = self.method(self.model_part, self.variable, self.container_location)

        list_of_values = GetListOfValues(values)
        list_of_indices = GetListOfValues(indices)

        list_of_str_values: 'list[str]' = []
        for i, v in enumerate(list_of_values):
            v_str = GetLengthAdjustedValue(v, self.header_lengths[i*2], self.precision)
            list_of_str_values.append(v_str)

            i_str = GetLengthAdjustedValue(list_of_indices[i], self.header_lengths[i*2+1], self.precision)
            list_of_str_values.append(i_str)

        return ", ".join(list_of_str_values)

    def GetHeadersString(self, header_details: 'dict[str, bool]') -> 'str':
        value_length = GetFormattedFloat(1e+100, self.precision)
        id_length = GetFormattedInt(int(1e+100), self.precision)

        list_of_headers = []
        for header in self._GetSuffixedTaggedHeaders(header_details, self.GetMethodInfo()):
            value_header = header.replace("<TAG>", "value")
            if len(value_header) < len(value_length):
                value_header = ("{:>" + str(len(value_length)) + "s}").format(value_header)
            list_of_headers.append(value_header)
            self.header_lengths.append(len(value_header))

            id_header = header.replace("<TAG>", "id")
            if len(id_header) < len(id_length):
                id_header = ("{:>" + str(len(id_length)) + "s}").format(id_header)
            list_of_headers.append(id_header)
            self.header_lengths.append(len(id_header))

        return ", ".join(list_of_headers)

class SpatialStatisticsVarianceOperation(SpatialStatisticsOperation):
    def __init__(self, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, precision: int) -> None:
        super().__init__(model_part, variable, container_location, norm)
        self.precision = precision
        self.header_lengths: 'list[int]' = []

    def GetMethodInfo(self) -> str:
        return "variance"

    def GetValueString(self) -> str:
        if self.norm is not None:
            means, variances = KratosStat.SpatialMethods.Variance(self.model_part, self.variable, self.container_location, self.norm)
        else:
            means, variances = KratosStat.SpatialMethods.Variance(self.model_part, self.variable, self.container_location)

        list_of_means = GetListOfValues(means)
        list_of_variances = GetListOfValues(variances)

        list_of_str_values: 'list[str]' = []
        for i, v in enumerate(list_of_means):
            v_mean_str = GetLengthAdjustedValue(v, self.header_lengths[i*2], self.precision)
            list_of_str_values.append(v_mean_str)

            v_variance_str = GetLengthAdjustedValue(list_of_variances[i], self.header_lengths[i*2+1], self.precision)
            list_of_str_values.append(v_variance_str)

        return ", ".join(list_of_str_values)

    def GetHeadersString(self, header_details: 'dict[str, bool]') -> 'str':
        value_length = GetFormattedFloat(1e+100, self.precision)

        list_of_headers = []
        for header in self._GetSuffixedTaggedHeaders(header_details, ""):
            mean_header = header.replace("<TAG>", "mean")
            if len(mean_header) < len(value_length):
                mean_header = ("{:>" + str(len(value_length)) + "s}").format(mean_header)
            list_of_headers.append(mean_header)
            self.header_lengths.append(len(mean_header))

            variance_header = header.replace("<TAG>", "variance")
            if len(variance_header) < len(value_length):
                variance_header = ("{:>" + str(len(value_length)) + "s}").format(variance_header)
            list_of_headers.append(variance_header)
            self.header_lengths.append(len(variance_header))

        return ", ".join(list_of_headers)

class SpatialStatisticsDistributionOperation(SpatialStatisticsOperation):
    def __init__(self, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, precision: int, parameters: Kratos.Parameters) -> None:
        super().__init__(model_part, variable, container_location, norm)
        self.precision = precision
        self.header_lengths: 'list[int]' = []
        self.parameters = parameters

        self.number_of_value_groups = 10
        if self.parameters.Has("number_of_value_groups"):
            self.number_of_value_groups = self.parameters["number_of_value_groups"].GetInt()

        if self.parameters.Has("min_value"):
            self.min_value = self.parameters["min_value"].GetDouble()
        else:
            raise RuntimeError("Using \"min\" for \"min_value\" is prohibited in this process.")

        if self.parameters.Has("max_value"):
            self.max_value = self.parameters["max_value"].GetDouble()
        else:
            raise RuntimeError("Using \"max\" for \"max_value\" is prohibited in this process.")

        self.group_upper_limits = []
        for i in range(self.number_of_value_groups + 1):
            self.group_upper_limits.append(self.min_value + (self.max_value - self.min_value) * i / self.number_of_value_groups)
        self.group_upper_limits.append(self.max_value)

    def GetMethodInfo(self) -> str:
        return "distribution"

    def GetValueString(self) -> str:
        if self.norm is not None:
            distribution_info = KratosStat.SpatialMethods.Distribution(self.model_part, self.variable, self.container_location, self.parameters, self.norm)
        else:
            distribution_info = KratosStat.SpatialMethods.Distribution(self.model_part, self.variable, self.container_location, self.parameters)

        min_value = distribution_info.GetMin()
        max_value = distribution_info.GetMax()
        group_number_of_values = distribution_info.GetGroupNumberOfValues()
        group_distribution = distribution_info.GetGroupValueDistributionPercentage()
        group_means = distribution_info.GetGroupMeans()
        group_variances = distribution_info.GetGroupVariances()

        list_of_min_values = GetListOfValues(min_value)
        list_of_max_values = GetListOfValues(max_value)

        group_details = []
        for i, _ in enumerate(self.group_upper_limits):
            list_of_number_of_values = GetListOfValues(group_number_of_values[i])
            list_of_distribution_values = GetListOfValues(group_distribution[i])
            for j, _ in enumerate(list_of_distribution_values):
                list_of_distribution_values[j] *= 100.0
            list_of_mean_values = GetListOfValues(group_means[i])
            list_of_variance_values = GetListOfValues(group_variances[i])
            group_details.append([list_of_number_of_values, list_of_distribution_values, list_of_mean_values, list_of_variance_values])

        values_per_component = 1 + 1 + len(self.group_upper_limits) * 4

        list_of_str_values: 'list[str]' = []
        for i, v_min in enumerate(list_of_min_values):
            v_min_str = GetLengthAdjustedValue(v_min, self.header_lengths[i*values_per_component], self.precision)
            list_of_str_values.append(v_min_str)

            v_max_str = GetLengthAdjustedValue(list_of_max_values[i], self.header_lengths[i*values_per_component+1], self.precision)
            list_of_str_values.append(v_max_str)

            # now add the group details
            index = i*values_per_component+2
            for group_detail in group_details:
                for sub_group_detail in group_detail:
                    v_str = GetLengthAdjustedValue(sub_group_detail[i], self.header_lengths[index], self.precision)
                    list_of_str_values.append(v_str)
                    index += 1

        return ", ".join(list_of_str_values)

    def GetHeadersString(self, header_details: 'dict[str, bool]') -> 'str':
        value_length = GetFormattedFloat(1e+100, self.precision)
        id_length = GetFormattedInt(int(1e+100), self.precision)

        list_of_headers = []
        for header in self._GetSuffixedTaggedHeaders(header_details, ""):
            min_header = header.replace("<TAG>", "min")
            if len(min_header) < len(value_length):
                min_header = ("{:>" + str(len(value_length)) + "s}").format(min_header)
            list_of_headers.append(min_header)
            self.header_lengths.append(len(min_header))

            max_header = header.replace("<TAG>", "max")
            if len(max_header) < len(value_length):
                max_header = ("{:>" + str(len(value_length)) + "s}").format(max_header)
            list_of_headers.append(max_header)
            self.header_lengths.append(len(max_header))

            # now add group details
            symbols = ["<" for _ in range(self.number_of_value_groups + 2)]
            symbols[-2] = "<="
            symbols[-1] = ">"
            for i, v in enumerate(self.group_upper_limits):
                v_str = GetFormattedFloat(v, self.precision)

                number_header = header.replace("<TAG>", f"grp_{i+1} # {symbols[i]} {v_str}")
                if len(number_header) < len(id_length):
                    number_header = ("{:>" + str(len(id_length)) + "s}").format(number_header)
                list_of_headers.append(number_header)
                self.header_lengths.append(len(number_header))

                distribution_header = header.replace("<TAG>", f"grp_{i+1} [%] {symbols[i]} {v_str}")
                if len(distribution_header) < len(value_length):
                    distribution_header = ("{:>" + str(len(value_length)) + "s}").format(distribution_header)
                list_of_headers.append(distribution_header)
                self.header_lengths.append(len(distribution_header))

                mean_header = header.replace("<TAG>", f"grp_{i+1} mean {symbols[i]} {v_str}")
                if len(mean_header) < len(value_length):
                    mean_header = ("{:>" + str(len(value_length)) + "s}").format(mean_header)
                list_of_headers.append(mean_header)
                self.header_lengths.append(len(mean_header))

                variance = header.replace("<TAG>", f"grp_{i+1} variance {symbols[i]} {v_str}")
                if len(variance) < len(value_length):
                    variance = ("{:>" + str(len(value_length)) + "s}").format(variance)
                list_of_headers.append(variance)
                self.header_lengths.append(len(variance))

        return ", ".join(list_of_headers)

def GetSpatialStatisticsOperation(method_name: str, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, precision: int, parameters: Kratos.Parameters) -> SpatialStatisticsOperation:
    method_dict = {
        "sum"           : lambda : SpatialStatisticsValueOperation(KratosStat.SpatialMethods.Sum, model_part, variable, container_location, norm, precision),
        "mean"          : lambda : SpatialStatisticsValueOperation(KratosStat.SpatialMethods.Mean, model_part, variable, container_location, norm, precision),
        "rootmeansquare": lambda : SpatialStatisticsValueOperation(KratosStat.SpatialMethods.RootMeanSquare, model_part, variable, container_location, norm, precision),
        "variance"      : lambda : SpatialStatisticsVarianceOperation(model_part, variable, container_location, norm, precision),
        "min"           : lambda : SpatialStatisticsValueIndexPairOperation(KratosStat.SpatialMethods.Min, model_part, variable, container_location, norm, precision),
        "max"           : lambda : SpatialStatisticsValueIndexPairOperation(KratosStat.SpatialMethods.Max, model_part, variable, container_location, norm, precision),
        "median"        : lambda : SpatialStatisticsValueIndexPairOperation(KratosStat.SpatialMethods.Median, model_part, variable, container_location, norm, precision),
        "distribution"  : lambda : SpatialStatisticsDistributionOperation(model_part, variable, container_location, norm, precision, parameters),
    }

    if method_name in method_dict.keys():
        return method_dict[method_name]()
    else:
        raise RuntimeError(f"Unsupported method name [ provided method name = \"{method_name}\" ]. Followings are supported:\n\t" + "\n\t".join(method_dict.keys()))
