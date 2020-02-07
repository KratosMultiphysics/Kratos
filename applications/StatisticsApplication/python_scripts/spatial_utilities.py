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

def GetNormTypeContainer(item_container, norm_type):
    if (norm_type == "value"):
        return item_container.ValueMethods
    else:
        return item_container.NormMethods

def GetMethod(item_norm_container, method_name):
    method_info_list = []
    for method in dir(item_norm_container):
        if (not method.startswith("__")):
            method_info_list.append([method.lower(), method])

    method_names_list = [method_info_list[i][0] for i in range(len(method_info_list))]
    method_list = [method_info_list[i][1] for i in range(len(method_info_list))]

    if (method_name not in method_names_list):
        msg = "Unknown method name [ \"method_name\" = \"" + method_name + "\" ]\n"
        msg += "Allowed method names are:\n    "
        msg += "\n    ".join(sorted(method_names_list))
        raise Exception(msg)

    return getattr(item_norm_container, method_list[method_names_list.index(method_name)])