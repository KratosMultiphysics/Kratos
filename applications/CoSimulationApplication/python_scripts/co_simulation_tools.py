import warnings

def ValidateAndAssignInputParameters(default, input, warnUnused=True):
    output = dict()
    for key in default.keys():
        if key in input.keys():

            value = input[key]
            #value = value.encode('utf-8')

            if type(default[key]) == type(value) or default[key] == type(value):
                output[key] = value
            else:
                error =  "ERROR: Key ", key, "has not the same type as the default value! "
                input_value = "INPUT: ", type(value)
                default_value = "DEFAULT: " + str(type(default[key])) + str(type(value))
                raise ValueError(error, input_value, default_value)

        else:
            if type(default[key]) == type:
                error =  "ERROR: mandatory key ", key, ":", default[key], "is missing! "
                raise ValueError(error)

            else:
                output[key] = default[key]
                warning_msg = "DEFAULT setting is used: " + str(key) + ":" + str(default[key])
                warnings.warn(warning_msg, Warning)

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

def GetSolvers(SolversDataList):
    solvers_map = {}
    num_solvers = len(SolversDataList)

    for i in range(0,num_solvers):
        solvers_module_name = "custom_co_simulation_solver_interfaces"
        __import__(solvers_module_name)
        import co_simulation_solver_factory as factory
        solver = factory.CreateSolverInterface(SolversDataList[i])
        solver_name = SolversDataList[i]["name"]
        solvers_map[solver_name] = solver

    return solvers_map


def GetSolverCoSimulationDetails(co_simulation_solver_settings):
    num_solvers = len(co_simulation_solver_settings)
    solver_cosim_details = {}
    for solver_settings in co_simulation_solver_settings:
        solver_name = solver_settings["name"]
        solver_cosim_details[solver_name] = solver_settings
    # TODO: check if the data is consitently defined! => maybe do at another place though...
    # - input in one is output in another
    # - one IO is defined for each data_name
    # - if the same data is defined multiple times
    # - check if data format has been specified
    return solver_cosim_details


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