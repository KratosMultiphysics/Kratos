from abc import ABC, abstractmethod
from typing import Any, Union
from math import log10

import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStat

def GetNorm(norm_name: str) -> Any:
    if norm_name == "noen":
        return None
    elif norm_name == "l2":
        return KratosStat.Norms.L2()
    elif norm_name == "infinity":
        return KratosStat.Norms.Infinity()
    elif norm_name == "trace":
        return KratosStat.Norms.Trace()
    elif norm_name.startswith("pnorm_"):
        return KratosStat.Norms.P(float(norm_name[6:]))
    elif norm_name.startswith("lpqnorm_("):
        pos = norm_name.rfind(",")
        if pos != -1:
            p = float(norm_name[9:pos])
            q = float(norm_name[pos+1:-1])
        else:
            raise RuntimeError(f"Unsupported norm info provided for LPQ norm [ provided info = \"{norm_name}\", required info = \"lpqnorm_(p,q)\" ].")
        return KratosStat.Norms.LPQ(p, q)
    else:
        raise RuntimeError(f"Unsupported norm name requested [ requested norm name = \"{norm_name}\" ]. Followings are supported:\n\tnone\n\tl2\n\tinfinity\n\ttrace\n\tpnorm_p\n\tlpqnorm_(p,q)")

def GetContainerLocation(container_location_name: str) -> Kratos.Globals.DataLocation:
    locations_dict = {
        "nodal_historical"        : Kratos.Globals.DataLocation.NodeHistorical,
        "nodal_non_historical"    : Kratos.Globals.DataLocation.NodeNonHistorical,
        "condition_non_historical": Kratos.Globals.DataLocation.Condition,
        "element_non_historical"  : Kratos.Globals.DataLocation.Element
    }

    if container_location_name in locations_dict.keys():
        return locations_dict[container_location_name]
    else:
        raise RuntimeError(f"Unsupported container location name [ provided container location name = \"{container_location_name}\" ]. Followings are supported:\n\t" + "\n\t".join(locations_dict.keys()))

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

def GetFormattedInt(value: int, precision: int) -> str:
    abs_value = abs(value)
    if abs_value >= 1e+100:
        exponents = int(log10(log10(abs_value))) + 1
    else:
        exponents = 2

    fixed_scientific_str = ("{: ." + str(precision - exponents) + "e}").format(value)
    fixed_str = ("{:> " + str(len(fixed_scientific_str)) + "d}").format(value)
    if len(fixed_str) <= len(fixed_scientific_str):
        return fixed_str
    else:
        return fixed_scientific_str

def GetFormattedFloat(value: float, precision: int) -> str:
    abs_value = abs(value)
    if abs_value >= 1e+100:
        exponents = int(log10(log10(abs_value))) + 1
    else:
        exponents = 2

    fixed_scientific_str = ("{: ." + str(precision - exponents) + "e}").format(value)
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
    @abstractmethod
    def GetValueString(self) -> str:
        pass

    @abstractmethod
    def GetHeadersString(self, header_details: 'dict[str, bool]') -> 'str':
        pass

class SpatialStatisticsValueOperation(SpatialStatisticsOperation):
    def __init__(self, method: Any, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, precision: int) -> None:
        self.model_part = model_part
        self.method = method
        self.variable = variable
        self.container_location = container_location
        self.norm = norm
        self.precision = precision
        self.header_lengths: 'list[int]' = []

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
        if self.norm is None:
            dummy_value = GetShapeSynchronizedValue(self.model_part, self.variable, self.container_location)
        else:
            dummy_value = 0.0

        value_length = GetFormattedFloat(1e+100, self.precision)

        list_of_headers = []
        for suffix in GetListOfValueHeaderSuffixes(dummy_value):
            header = GetTaggedHeader(header_details, self.container_location.name, self.variable.Name(), self.norm, suffix, self.method.__name__.lower())
            header = header.replace("<TAG>", "").strip()
            if len(header) < len(value_length):
                header = ("{:>" + str(len(value_length)) + "s}").format(header)
            list_of_headers.append(header)
            self.header_lengths.append(len(header))

        return ", ".join(list_of_headers)

class SpatialStatisticsValueIndexPairOperation(SpatialStatisticsOperation):
    def __init__(self, method: Any, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, precision: int) -> None:
        self.model_part = model_part
        self.method = method
        self.variable = variable
        self.container_location = container_location
        self.norm = norm
        self.precision = precision
        self.header_lengths: 'list[int]' = []

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
        if self.norm is None:
            dummy_value = GetShapeSynchronizedValue(self.model_part, self.variable, self.container_location)
        else:
            dummy_value = 0.0

        value_length = GetFormattedFloat(1e+100, self.precision)
        id_length = GetFormattedInt(int(1e+100), self.precision)

        list_of_headers = []
        for suffix in GetListOfValueHeaderSuffixes(dummy_value):
            header = GetTaggedHeader(header_details, self.container_location.name, self.variable.Name(), self.norm, suffix, self.method.__name__.lower())

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
        self.model_part = model_part
        self.variable = variable
        self.container_location = container_location
        self.norm = norm
        self.precision = precision
        self.header_lengths: 'list[int]' = []

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
        if self.norm is None:
            dummy_value = GetShapeSynchronizedValue(self.model_part, self.variable, self.container_location)
        else:
            dummy_value = 0.0

        value_length = GetFormattedFloat(1e+100, self.precision)

        list_of_headers = []
        for suffix in GetListOfValueHeaderSuffixes(dummy_value):
            header = GetTaggedHeader(header_details, self.container_location.name, self.variable.Name(), self.norm, suffix, "")

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
        self.model_part = model_part
        self.variable = variable
        self.container_location = container_location
        self.norm = norm
        self.precision = precision
        self.parameters = parameters
        self.header_lengths: 'list[int]' = []

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

    def GetValueString(self) -> str:
        if self.norm is not None:
            distribution_info: KratosStat.SpatialMethods.Array3DistributionInfo = KratosStat.SpatialMethods.Distribution(self.model_part, self.variable, self.container_location, self.parameters, self.norm)
        else:
            distribution_info: KratosStat.SpatialMethods.Array3DistributionInfo = KratosStat.SpatialMethods.Distribution(self.model_part, self.variable, self.container_location, self.parameters)

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
        if self.norm is None:
            dummy_value = GetShapeSynchronizedValue(self.model_part, self.variable, self.container_location)
        else:
            dummy_value = 0.0

        value_length = GetFormattedFloat(1e+100, self.precision)
        id_length = GetFormattedInt(int(1e+100), self.precision)

        list_of_headers = []
        for suffix in GetListOfValueHeaderSuffixes(dummy_value):
            header = GetTaggedHeader(header_details, self.container_location.name, self.variable.Name(), self.norm, suffix, "")

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

if __name__ == "__main__":
    model = Kratos.Model()
    model_part = model.CreateModelPart("test")
    model_part.CreateNewNode(1, 1, 1, 1)
    model_part.CreateNewNode(2, 2, 2, 2)
    for node in model_part.Nodes:
        node.SetValue(Kratos.VELOCITY, Kratos.Array3([-node.Id*1e90, node.Id, node.Id*0]))

    print("---------------------------------------------")
    output = GetSpatialStatisticsOperation("sum", model_part, Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeNonHistorical, KratosStat.Norms.L2(), 12, None)
    print(output.GetHeadersString({"location": True, "norm": False, "variable": False, "method": False}))
    print(output.GetValueString())
    print(output.GetValueString())
    print(output.GetValueString())

    print("----------------------Min-----------------------")
    output = GetSpatialStatisticsOperation("min", model_part, Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeNonHistorical, None, 12, None)
    print(output.GetHeadersString({"location": True, "norm": False, "variable": False, "method": False}))
    print(output.GetValueString())
    print(output.GetValueString())
    print(output.GetValueString())

    print("---------------------------------------------")
    output = GetSpatialStatisticsOperation("variance", model_part, Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeNonHistorical, None, 12, None)
    print(output.GetHeadersString({"location": True, "norm": False, "variable": False, "method": True}))
    print(output.GetValueString())
    print(output.GetValueString())
    print(output.GetValueString())

    print("---------------------------------------------")
    parameters = Kratos.Parameters("""{
        "min_value": -1.0,
        "max_value": 10.0
    }""")
    output = GetSpatialStatisticsOperation("distribution", model_part, Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeNonHistorical, KratosStat.Norms.L2(), 12, parameters)
    print(output.GetHeadersString({"location": True, "norm": False, "variable": False, "method": True}))
    print(output.GetValueString())
    print(output.GetValueString())
    print(output.GetValueString())

#         header += f" @ {self.container_location.name} {method_name}\""


#         if self.norm is not None:
#             return [method_name]
#         else:
#             container_retrieval_method, entity_data_retireval_method = GetDataRetrievalMethods(self.container_location)
#             v = GetValueForVariable(self.variable)
#             for entity in container_retrieval_method(self.model_part):
#                 v = entity_data_retireval_method(entity, self.variable)
#                 break
#             self.model_part.GetCommunicator().GetDataCommunicator().SynchronizeShape(v)
#             return [f"{method_name}{suffix}" for suffix in GetListOfValueHeaderSuffixes(v)]

# class SpatialStatisticsValueIndexPairOperation(SpatialStatisticsOperation):
#     def __init__(self, method: Any, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any) -> None:
#         self.model_part = model_part
#         self.method = method
#         self.variable = variable
#         self.container_location = container_location
#         self.norm = norm

#     def GetValues(self) -> 'list[Union[int, float]]':
#         if self.norm is not None:
#             data_value, data_index = self.method(self.model_part, self.variable, self.container_location, self.norm)
#             return [data_value, data_index]
#         else:
#             data_values, data_indices = self.method(self.model_part, self.variable, self.container_location)
#             values: 'list[Union[int, float]]' = []
#             for data_v, data_i in zip(GetListOfValues(data_values), GetListOfValues(data_indices)):
#                 values.append(data_v)
#                 values.append(data_i)
#             return values

#     def GetHeaders(self) -> 'list[str]':
#         method_name = self.method.__name__.lower()
#         if self.norm is not None:
#             return [f"{method_name}_value", f"{method_name}_id"]
#         else:
#             container_retrieval_method, entity_data_retireval_method = GetDataRetrievalMethods(self.container_location)
#             v = GetValueForVariable(self.variable)
#             for entity in container_retrieval_method(self.model_part):
#                 v = entity_data_retireval_method(entity, self.variable)
#                 break
#             self.model_part.GetCommunicator().GetDataCommunicator().SynchronizeShape(v)
#             headers: 'list[str]' = []
#             for suffix in GetListOfValueHeaderSuffixes(v):
#                 headers.append(f"{method_name}_value{suffix}")
#                 headers.append(f"{method_name}_id{suffix}")
#             return headers

# class SpatialStatisticsVarianceOperation(SpatialStatisticsOperation):
#     def __init__(self, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any) -> None:
#         self.model_part = model_part
#         self.variable = variable
#         self.container_location = container_location
#         self.norm = norm

#     def GetValues(self) -> 'list[Union[int, float]]':
#         if self.norm is not None:
#             mean, variance = KratosStat.SpatialMethods.Variance(self.model_part, self.variable, self.container_location, self.norm)
#             return [mean, variance]
#         else:
#             means, variances = KratosStat.SpatialMethods.Variance(self.model_part, self.variable, self.container_location)
#             values: 'list[Union[int, float]]' = []
#             for mean, variance in zip(GetListOfValues(means), GetListOfValues(variances)):
#                 values.append(mean)
#                 values.append(variance)
#             return values

#     def GetHeaders(self) -> 'list[str]':
#         if self.norm is not None:
#             return ["mean", "variance"]
#         else:
#             container_retrieval_method, entity_data_retireval_method = GetDataRetrievalMethods(self.container_location)
#             v = GetValueForVariable(self.variable)
#             for entity in container_retrieval_method(self.model_part):
#                 v = entity_data_retireval_method(entity, self.variable)
#                 break
#             self.model_part.GetCommunicator().GetDataCommunicator().SynchronizeShape(v)
#             headers: 'list[str]' = []
#             for suffix in GetListOfValueHeaderSuffixes(v):
#                 headers.append(f"mean{suffix}")
#                 headers.append(f"variance{suffix}")
#             return headers

# class SpatialStatisticsDistributionOperation(SpatialStatisticsOperation):
#     def __init__(self, model_part: Kratos.ModelPart, variable: Any, container_location: Kratos.Globals.DataLocation, norm: Any, parameters: Kratos.Parameters) -> None:
#         self.model_part = model_part
#         self.variable = variable
#         self.container_location = container_location
#         self.norm = norm
#         self.parameters = parameters

#     def GetValues(self) -> 'list[Union[int, float]]':
#         if self.norm is not None:
#             distribution_data = KratosStat.SpatialMethods.Distribution(self.model_part, self.variable, self.container_location, self.parameters, self.norm)
#             return [mean, variance]
#         else:
#             distribution_data = KratosStat.SpatialMethods.Distribution(self.model_part, self.variable, self.container_location, self.parameters)
#             values: 'list[Union[int, float]]' = []
#             for mean, variance in zip(GetListOfValues(means), GetListOfValues(variances)):
#                 values.append(mean)
#                 values.append(variance)
#             return values

#     def GetHeaders(self) -> 'list[str]':
#         if self.norm is not None:
#             return ["mean", "variance"]
#         else:
#             container_retrieval_method, entity_data_retireval_method = GetDataRetrievalMethods(self.container_location)
#             v = GetValueForVariable(self.variable)
#             for entity in container_retrieval_method(self.model_part):
#                 v = entity_data_retireval_method(entity, self.variable)
#                 break
#             self.model_part.GetCommunicator().GetDataCommunicator().SynchronizeShape(v)
#             headers: 'list[str]' = []
#             for suffix in GetListOfValueHeaderSuffixes(v):
#                 headers.append(f"mean{suffix}")
#                 headers.append(f"variance{suffix}")
#             return headers

# class SpatialStatisticsMeanOperation(SpatialStatisticsOperation):
#     def _GetMethod(self) -> Any:
#         return lambda *args: KratosStat.SpatialMethods.Mean(self.model_part, self.variable, self.container_location, *args)

# class SpatialStatisticsRootMeanSquareOperation(SpatialStatisticsOperation):
#     def _GetMethod(self) -> Any:
#         return lambda *args: KratosStat.SpatialMethods.RootMeanSquare(self.model_part, self.variable, self.container_location, *args)

# class SpatialStatisticsVarianceOperation(SpatialStatisticsOperation):
#     def _GetMethod(self) -> Any:
#         return lambda *args: KratosStat.SpatialMethods.Variance(self.model_part, self.variable, self.container_location, *args)

# class SpatialStatisticsMinOperation(SpatialStatisticsOperation):
#     pass

# class SpatialStatisticsMaxOperation(SpatialStatisticsOperation):
#     pass

# class SpatialStatisticsMedianOperation(SpatialStatisticsOperation):
#     pass

# class SpatialStatisticsDistributionOperation(SpatialStatisticsOperation):
#     pass



#     if method_name in method_dict.keys():
#         return method_dict[method_name]
#     else:
#         raise RuntimeError(f"Unsupported container location name [ provided container location name = \"{method_name}\" ]. Followings are supported:\n\t" + "\n\t".join(method_dict.keys()))


# class SpatialStatisticsProcess(Kratos.OutputProcess):
#     """A process to calculate spatial statistics on Kratos containers

#     This process calculates spatial statistics for given variables in a given container.

#     This process is compatible with OpenMP and MPI with restart

#     Args:
#         model (Kratos.Model): Model used in problem
#         settings (Kratos.Parameters): Kratos parameter settings for process
#     """
#     def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
#         Kratos.OutputProcess.__init__(self)

#         default_parameters = Kratos.Parameters("""
#         {
#             "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
#             "echo_level"              : 0,
#             "computation_processes"   : [],
#             "computation_points"      : [
#                 "ExecuteInitialize",
#                 "ExecuteInitializeSolutionStep",
#                 "Execute",
#                 "ExecuteFinalizeSolutionStep"
#             ],
#             "input_variable_settings" : [
#                 {
#                     "variable_names" : [],
#                     "norm_type"      : "none",
#                     "container"      : "nodal_historical"
#                 }
#             ],
#             "statistics_methods": [
#                 {
#                     "method_name"    : "",
#                     "method_settings": {}
#                 }
#             ],
#             "output_settings" : {
#                 "interval"               : [0.0, "End"],
#                 "output_control_variable": "STEP",
#                 "output_time_interval"   : 1,
#                 "write_kratos_version"   : false,
#                 "write_time_stamp"       : false,
#                 "output_value_precision" : 5,
#                 "output_value_length"    : 14,
#                 "output_file_settings"   : {
#                     "file_name"  : "<model_part_name>.dat",
#                     "output_path": "spatial_statistics_output",
#                     "write_buffer_size" : -1
#                 }
#             }
#         }  """)

#         self.model = model
#         self.settings = settings
#         self.settings.ValidateAndAssignDefaults(default_parameters)

#         self.operations = []

