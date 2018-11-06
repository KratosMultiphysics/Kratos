import warnings
import co_simulation_data_structure
cs_data_structure = co_simulation_data_structure.__KRATOS_DATA_STRUCTURE__
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

def ImportDataStructure(parameters_file_name):
    import json
    import sys
    global cs_data_structure
    with open(parameters_file_name,'r') as parameter_file:
        parameters = json.load(parameter_file) # This is for getting the flag for database
        if 'data_structure' in parameters['problem_data'].keys():
            data_structure = parameters['problem_data']['data_structure']
        else:
            data_structure = "co_sim_app"

        co_simulation_data_structure.Initialize(data_structure)
        cs_data_structure = co_simulation_data_structure.__KRATOS_DATA_STRUCTURE__

    return cs_data_structure
