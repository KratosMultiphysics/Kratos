def GetNormTypeContainer(item_container, norm_type):
    if (norm_type == "none" or norm_type == "value"):
        return item_container.ValueMethods
    else:
        return item_container.NormMethods

def GetAvailableMethods(item_norm_container):
    method_info_list = []
    for method in dir(item_norm_container):
        if (not method.startswith("__")):
            method_info_list.append([method.lower(), method])

    method_names_list = [method_info[0] for method_info in method_info_list]
    method_list = [method_info[1] for method_info in method_info_list]

    return method_names_list, method_list

def GetMethod(item_norm_container, method_name):
    method_names_list, method_list = GetAvailableMethods(item_norm_container)

    if method_name not in method_names_list:
        msg = "Unknown method name requested for the container " + str(item_norm_container) + ". [ \"method_name\" = \"" + method_name + "\" ]\n"
        msg += "Allowed method names are:\n    "
        msg += "\n    ".join(sorted(method_names_list))
        raise Exception(msg)

    return getattr(item_norm_container,
                   method_list[method_names_list.index(method_name)])
