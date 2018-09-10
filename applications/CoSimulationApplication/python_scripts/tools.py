import warnings

def SolverFactory(SolversParameters):
    SolversList = {}
    num_solvers = SolversDataList.size()

    for i in range(0,num_solvers):
        solver_module = __import__(SolversDataList[i]['type'].GetString())
        solver_name = SolversDataList[i]['name'].GetString()
        solver = solver_module.Create(SolversDataList[i])
        SolversList[solver_name] = solver

    return SolversList

def ConvergenceAcceleratorFactory(SolversDataList):
    SolversList = {}
    num_solvers = SolversDataList.size()

    for i in range(0,num_solvers):
        solver_module = __import__(SolversDataList[i]['type'].GetString())
        solver_name = SolversDataList[i]['name'].GetString()
        solver = solver_module.Create(SolversDataList[i])
        SolversList[solver_name] = solver

    return SolversList

def ConvergenceCriteriaFactory(SolversDataList):
    SolversList = {}
    num_solvers = SolversDataList.size()

    for i in range(0,num_solvers):
        solver_module = __import__(SolversDataList[i]['type'].GetString())
        solver_name = SolversDataList[i]['name'].GetString()
        solver = solver_module.Create(SolversDataList[i])
        SolversList[solver_name] = solver

    return SolversList    

def ValidateAndAssignInputParameters(default, input, warnUnused=True):
    output = dict()
    for key in default.keys():
        if key in input.keys():

            value = input[key]

            if type(value) in [unicode]:
                value = EncodeString(value)

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
                    value = EncodeString(value)
                output[inputSetting] = value

                if warnUnused:
                    warning_msg = "UNUSED setting in input: " + str(inputSetting) + ":" + str(value)
                    warnings.warn(warning_msg, Warning)

    return output