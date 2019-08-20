import sys


cs_data_structure = None


# ImportDataStructure: Imports the data structure which is specified in the parameters file
#
#  @param parameters_file_name   The JSON file name which contains the settings for the co-simulation
def ImportDataStructure(parameters_file_name):
    import json

    global cs_data_structure

    # Set default data_structure
    data_structure_name = 'pyKratos'
    # Read data_structure from parameter file if present
    with open(parameters_file_name, 'r') as parameter_file:
        parameters = json.load(parameter_file)
        if 'data_structure' in parameters['settings'].keys():
            data_structure_name = parameters['settings']['data_structure']

    # Initialize cs_data_structure and import corresponding module
    if cs_data_structure is None:
        if data_structure_name == 'KratosMultiphysics':
            cs_data_structure = __import__('KratosMultiphysics')
        elif data_structure_name == 'pyKratos':
            sys.path.append("../custom_data_structure/")
            cs_data_structure = __import__('pyKratos')
        else:
            raise Exception('data_structure needs to be KratosMultiphysics or pyKratos')

    return cs_data_structure


# CreateInstance: Creates an instance of a given class based on a settings dictionary
#
#  @param settings  The settings dictionary, with the key "type"
def CreateInstance(settings):
    object_type = settings["type"].GetString()
    object_module_full = "KratosMultiphysics.CoSimulationApplication."+object_type
    object_module = __import__(object_module_full, fromlist=[object_type])
    return object_module.Create(settings)


# InnerProduct: Computes the inner product for two given vectors (as python lists)
#
#  @param list_one   First vector (as a list)
#  @param list_two   Second vector (as a list)
def InnerProduct(list_one, list_two):
    result = 0.0
    num_entries = len(list_one)
    if len(list_one) == len(list_two):
        for i in range(num_entries):
            component_one = list_one[i]
            component_two = list_two[i]
            result += component_one * component_two

    return result


# CalculateNorm: Calculates the L2-norm of the list
#
#  @param list   Vector (as a list)
def CalculateNorm(list):
    norm = 0.0
    for value in list:
        norm += value*value

    return norm


# Class contains definition of colors. This is to be used as a struct.
#
# Example usage: print(bcolors.HEADER + "This is a header in header color" + bcolor.ENDC)
# IMPORTANT: The end of the print statement should always contain bcolor.ENDC
class bcolors:
    HEADER    = '\033[95m'
    BLUE      = '\033[94m'
    GREEN     = '\033[92m'
    MEGENTA   = '\033[96m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'


# PrintInfo: Printing information with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintInfo(label, *args):
    print(label, " ".join(map(str, args)))


# PrintWarning: Printing a warning with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintWarning(label, *args):
    print(label, " ".join(map(str, args)))
