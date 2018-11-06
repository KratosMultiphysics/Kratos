import warnings
import co_simulation_data_structure

# TODO: make it recursive by default
def ValidateAndAssignInputParameters(default, input, warnUnused=True):
    output = dict()
    for key in default.keys():
        if key in input.keys():
            value = input[key]
            if type(default[key]) == type(value) or default[key] == type(value):
                output[key] = value
            else:
                error =  "ERROR: Key ", key, "has not the same type as the default value! "
                input_value = "INPUT: ", type(value)
                default_value = "DEFAULT: " + str(type(default[key])) + str(type(value))
                raise ValueError(error, input_value, default_value)

        else:
            if type(default[key]) == type:
                error =  "ERROR: mandatory setting ", key, ":", default[key], "is missing! "
                raise ValueError(error)
            else:
                output[key] = default[key]

        for inputSetting in input.keys():
            if inputSetting not in default.keys():
                value = input[inputSetting]
                if type(value) == str:
                    value = value.encode('utf-8')
                output[inputSetting] = value

                if warnUnused:
                    warning_msg = "UNUSED setting in input: " + str(inputSetting) + ":" + str(value)
                    warnings.warn(warning_msg, Warning)

    return output



## Class contains definition of colors. This is to be used as a struct
#
# Example usage print(bcolors.HEADER + "This is a header in header color" + bcolor.ENDC)
# IMPORTANT : The end of the print statement should always containt bcolor.ENDC
class bcolors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    MEGENTA = '\033[96m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def MakeDataConfig(custom_config):
    defaultSettings = {}
    defaultSettings["name"] = str # MANDATORY
    defaultSettings["format"] = "python_list"
    defaultSettings["dimension"] = int # MANDATORY
    defaultSettings["size"] = 0
    defaultSettings["mesh_name"] = ""
    defaultSettings["location_on_mesh"] = ""
    settings = ValidateAndAssignInputParameters(defaultSettings, custom_config, False)

    if settings["format"] == "python_list" :
        settings["data"] = []

    return settings


def ImportDataStructure(parameters_file_name):
    import json
    import sys
    with open(parameters_file_name,'r') as parameter_file:
        parameters = json.load(parameter_file) # This is for getting the flag for database
        if 'data_structure' in parameters['problem_data'].keys():
            data_structure = parameters['problem_data']['data_structure']
        else:
            data_structure = "co_sim_app"

        co_simulation_data_structure.Initialize(data_structure)
        cs_data_structure = co_simulation_data_structure.__KRATOS_DATA_STRUCTRURE__

    return cs_data_structure
