import math
cs_data_structure = None

## Class contains definition of colors. This is to be used as a struct
#
# Example usage print(bcolors.HEADER + "This is a header in header color" + bcolor.ENDC)
# IMPORTANT : The end of the print statement should always contain bcolor.ENDC
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


## ImportDataStructure : Imports the data structure which is specified in the parameters file
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
        if 'data_structure' in parameters['problem_data'].keys():
            data_structure_name = parameters['problem_data']['data_structure']

    # Initialize cs_data_structure and import corresponding module
    if cs_data_structure is None:
        if data_structure_name == 'KratosMultiphysics':
            cs_data_structure = __import__('KratosMultiphysics')
        elif data_structure_name == 'pyKratos':
            cs_data_structure = __import__('pyKratos')
        else:
            raise Exception('data_structure needs to be KratosMultiphysics or pyKratos')

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
    data_def = solver.data_map[data_name]
    data_mesh = solver.model[data_def["geometry_name"].GetString()]
    data_variable = cs_data_structure.KratosGlobals.GetVariable(data_name)
    for node in data_mesh.Nodes:
        data_value = node.GetSolutionStepValue(data_variable,0) #TODO what if non-historical?
        for value in data_value:
            data.append(value)

    return data

## ApplyUpdateToData : Apply the update provided to the data with name data_name
#
#  @param solver        The solver from which data is to be obtained.
#  @param data_name     The name of the data to which the update is to be applied
#  @param update        The update list which is to be applied to the data with name data_name
def ApplyUpdateToData(solver, data_name, updated_data):
    data_def = solver.data_map[data_name]
    data_mesh = solver.model[data_def["geometry_name"].GetString()]
    data_variable = cs_data_structure.KratosGlobals.GetVariable(data_name)
    index = 0
    for node in data_mesh.Nodes: # #TODO local nodes to also work in MPI?
        updated_value = []
        value = node.GetSolutionStepValue(data_variable,0)
        # TODO @aditya the data might also be non-historical => GetValue
        for i, value_i in enumerate(value):
            #updated_value.append(0.0)
            updated_value.append(updated_data[index])
            index = index + 1

        node.SetSolutionStepValue(data_variable, 0, updated_value)

## PrintInfo : Printing information with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintInfo(label, *args):
    print(label, " ".join(map(str,args)))

## PrintInfo : Printing a warning with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintWarning(label, *args):
    print(label, " ".join(map(str,args)))



## Class CouplingInterfaceData: Class to hold different properties of the data field contributed in
#                           CoSimulation.
#
class CouplingInterfaceData(object):
    def __init__(self, custom_config, solver):

        default_config = cs_data_structure.Parameters("""
        {
            "name" : "default",
            "dimension" : 0,
            "geometry_name" : "",
            "location_on_mesh":"on_nodes"
        }
        """)
        custom_config.ValidateAndAssignDefaults(default_config)

        self.name = custom_config["name"].GetString()
        self.variable = None
        self.filters = []
        self.solver = solver
        self.dimension = custom_config["dimension"].GetInt()
        self.location_on_mesh = custom_config["location_on_mesh"].GetString()
        self.mesh_name = custom_config["geometry_name"].GetString()
        self.destination_data = None
        self.origin_data = None
        self.mapper_settings = None

    def ApplyFilters(self):
        for filter in self.filters:
            #filter.InitializeNonLinearIteration()
            filter.Apply()
            #filter.FinalizeNonLinearIteration()

    def GetPythonList(self):
        data = []
        data_mesh = self.solver.model[self.mesh_name]
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        for node in data_mesh.Nodes:
            data_value = node.GetSolutionStepValue(data_variable,0) #TODO what if non-historical?
            for value in data_value:
                data.append(value)
        return data

    def GetNumpyArray(self):
        pass

    def ApplyUpdateToData(self, update):
        data_mesh = self.solver.model[self.mesh_name]
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        index = 0
        for node in data_mesh.Nodes: # #TODO: local nodes to also work in MPI?
            updated_value = []
            value = node.GetSolutionStepValue(data_variable,0)
            # TODO: aditya the data might also be non-historical => GetValue
            for value_i in value:
                updated_value.append(update[index])
                index = index + 1
            node.SetSolutionStepValue(data_variable, 0, updated_value)


