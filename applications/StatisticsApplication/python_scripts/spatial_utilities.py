import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as Statistics

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

def GetMethodHeaders(method_name, parameters):
    method_headers = [
        ["sum", ["_Sum"]],
        ["mean", ["_Mean"]],
        ["rootmeansquare", ["_RootMeanSquare"]],
        ["variance", ["_Mean", "_Variance"]],
        ["min", ["_Min", "_Id"]],
        ["max", ["_Max", "_Id"]],
        ["median", ["_Median"]],
        ["distribution", ["_Min", "_Max"]]
    ]

    if (method_name == "distribution"):
        method_headers_list = method_headers[6][1]
        number_of_groups = parameters["number_of_value_groups"].GetInt()
        method_headers_list.append(" < Min")
        for i in range(number_of_groups-1):
            method_headers_list.append(" < " + str((i+1) / number_of_groups))
        method_headers_list.append(" <= Max")
        method_headers_list.append(" > Max")

    method_header_names = [ method_headers[i][0] for i in range(len(method_headers)) ]
    return method_headers[method_header_names.index(method_name)][1]

def GetMethodValues(method_name, norm_type, variable_name, output):
    if norm_type == "none":
        variable_sub_list = GetVariableHeaders(norm_type, variable_name)
        if (method_name in ["sum", "mean", "rootmeansquare"]):
            if (len(variable_sub_list) == 1):
                return str(output)
            else:
                str_output = ""
                for i in range(len(variable_sub_list)):
                    str_output += str(output[i]) + " "
                str_output = str_output[:-1]
                return str_output
        elif (method_name in ["variance"]):
            if (len(variable_sub_list) == 1):
                return str(output[0]) + " " + str(output[1])
            else:
                str_output = ""
                for i in range(len(variable_sub_list)):
                    str_output += str(output[0][i]) + " "
                for i in range(len(variable_sub_list)):
                    str_output += str(output[1][i]) + " "
                str_output = str_output[:-1]
                return str_output
        else:
            raise Exception("Unknown value method " + method_name)
    else:
        if (method_name in ["sum", "mean", "median", "rootmeansquare"]):
            return str(output)
        elif (method_name in ["min", "max", "variance"]):
            return str(output[0]) + " " + str(output[1])
        elif (method_name == "distribution"):
            msg = str(output[0]) + " " + str(output[1])
            for _v in output[4]:
                msg += " " + str(_v)
            return msg
        else:
            raise Exception("Unknown norm method " + method_name)

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

