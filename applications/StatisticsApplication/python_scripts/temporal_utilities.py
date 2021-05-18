import KratosMultiphysics.StatisticsApplication as Statistics

def GetItemContainer(item_container_name):
    item_container_types = [
        [
            "nodal_historical_historical",
            Statistics.TemporalMethods.Historical.HistoricalOutput
        ],
        [
            "nodal_historical_non_historical",
            Statistics.TemporalMethods.Historical.NonHistoricalOutput
        ],
        [
            "nodal_non_historical",
            Statistics.TemporalMethods.NonHistorical.Nodes
        ],
        [
            "element_non_historical",
            Statistics.TemporalMethods.NonHistorical.Elements
        ],
        [
            "condition_non_historical",
            Statistics.TemporalMethods.NonHistorical.Conditions
        ],
    ]

    item_container_names_list = [
        item_container_types[i][0] for i in range(len(item_container_types))
    ]
    item_container_types_list = [
        item_container_types[i][1] for i in range(len(item_container_types))
    ]

    if (item_container_name not in item_container_names_list):
        msg = "Unknown container [ \"container\" = \"" + item_container_name + "\" ]\n"
        msg += "Allowed containers are:\n    "
        msg += "\n    ".join(sorted(item_container_names_list))
        raise Exception(msg)

    return item_container_types_list[item_container_names_list.index(
        item_container_name)]


def GetMethodClass(item_norm_method_container, variable_type):
    method_class_info_list = []
    for method in dir(item_norm_method_container):
        if (not method.startswith("__")):
            method_class_info_list.append([method.lower(), method])

    variable_types_list = [method_class_info_list[i][0] for i in range(len(method_class_info_list))]
    method_type_class_list = [method_class_info_list[i][1] for i in range(len(method_class_info_list))]

    if (variable_type not in variable_types_list):
        msg = "Unknown variable type. [ variable_type = " + variable_type + " ]\n"
        msg += "Allowed variable types are:\n    "
        msg += "\n    ".join(sorted(variable_types_list))
        raise Exception(msg)

    return getattr(item_norm_method_container, method_type_class_list[variable_types_list.index(variable_type)])



