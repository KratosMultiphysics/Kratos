import warnings
import co_simulation_data_structure
cs_data_structure = co_simulation_data_structure.__KRATOS_DATA_STRUCTURE__
import math
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


## ImportDataStructure : Imports the data structre which is specified in the parameters file
#
#  @param parameters_file_name   The JSON file name which contains the settings for the co-simulation
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

## InnerProduct : Computes the inner product for two give vectors (as python lists)
#
#  @param list_one   First vector (as a list)
#  @param list_two   Second vector (as a list)
def InnerProduct(list_one, list_two):
    result = 0.0
    num_entries = len(list_one)
    if(len(list_one) == len(list_two)):
        for i in range(num_entries):
            component_one = list_one[i]
            component_two = list_two[i]
            result += component_one * component_two

    return result

## CalculateNorm : Calculates the two norm of the list (residual) passed in
#                       Variants of this class can use numpy for calculating this norm.
#
#  @param self            The object pointer
#  @param residual_list   The list of the residual
def CalculateNorm(residual_list): ## Can use numpy here.
    norm = 0.0
    num_entries = len(residual_list)
    for residual in residual_list:
        norm += residual * residual

    norm = math.sqrt(norm)/num_entries
    return norm


## GetDataAsList : Converts the data as a python list.
#                       variants of this class can use numpy array
#
#  @param solver        The solver from which data is to be obtained.
#  @param data_name     The name of the data which is to be obtained as a list
def GetDataAsList(solver, data_name):
    data = []
    data_def = solver.data_list[data_name]
    data_mesh = solver.model[data_def["geometry_name"].GetString()]
    data_variable = cs_data_structure.KratosGlobals.GetVariable(data_name)
    for node in data_mesh.Nodes:
        data_value = node.GetSolutionStepValue(data_variable,0)
        for value in data_value:
            data.append(value)

    return data

## ApplyUpdateToData : Apply the update provided to the data with name data_name
#
#  @param solver        The solver from which data is to be obtained.
#  @param data_name     The name of the data to which the update is to be applied
#  @param update        The update list which is to be applied to the data with name data_name
def ApplyUpdateToData(solver, data_name, update):
    data_def = solver.data_list[data_name]
    data_mesh = solver.model[data_def["geometry_name"].GetString()]
    data_variable = cs_data_structure.KratosGlobals.GetVariable(data_name)
    index = 0

    for node in data_mesh.Nodes:
        updated_value = []
        value = node.GetSolutionStepValue(data_variable,0)
        for i, value_i in enumerate(value):
            updated_value.append(value_i + update[index])
        node.SetSolutionStepValue(data_variable, updated_value)
        index = index + 1