import KratosMultiphysics
import KratosMultiphysics.mpi

def GetVariableAndType(varName, splitByCommponents):
    """ This method returns the variable with the given name

    Keyword arguments:
    varName -- The name of the variable to return

    """

    KratosGlobals = KratosMultiphysics.KratosGlobals

    if KratosGlobals.Kernel.HasDoubleVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetDoubleVariable(varName)}
    elif KratosGlobals.Kernel.HasArrayVariable(varName):
        if splitByCommponents:
            for component in ["X","Y","Z","XX","YY","XY"]:
                if KratosGlobals.Kernel.HasDoubleVariable(varName+"_"+component):
                    yield {"name": varName+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(varName+"_"+component)}
        else:
            yield {"name": varName, "value": KratosGlobals.Kernel.GetArrayVariable(varName)}
    elif KratosGlobals.Kernel.HasArray4Variable(varName):
        if splitByCommponents:
            for component in ["XX","XY","YX","YY"]:
                yield {"name": varName+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(varName+"_"+component)}
        else:
            yield {"name": varName, "value": KratosGlobals.Kernel.GetArray4Variable(varName)}
    elif KratosGlobals.Kernel.HasArray6Variable(varName):
        if splitByCommponents:
            for component in ["XX","YY","ZZ","XY","YZ","XZ"]:
                yield {"name": varName+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(varName+"_"+component)}
        else:
            yield {"name": varName, "value": KratosGlobals.Kernel.GetArray6Variable(varName)}
    elif KratosGlobals.Kernel.HasArray9Variable(varName):
        if splitByCommponents:
            for component in ["XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ"]:
                yield {"name": varName+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(varName+"_"+component)}
        else:
            yield {"name": varName, "value": KratosGlobals.Kernel.GetArray9Variable(varName)}
    elif KratosGlobals.Kernel.HasBoolVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetBoolVariable(varName)}
    elif KratosGlobals.Kernel.HasIntVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetIntVariable(varName)}
    elif KratosGlobals.Kernel.HasUnsignedIntVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetUnsignedIntVariable(varName)}
    elif KratosGlobals.Kernel.HasVectorVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetVectorVariable(varName)}
    elif KratosGlobals.Kernel.HasMatrixVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetMatrixVariable(varName)}
    elif KratosGlobals.Kernel.HasStringVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetStringVariable(varName)}
    elif KratosGlobals.Kernel.HasFlagsVariable(varName):
        yield {"name": varName, "value": KratosGlobals.Kernel.GetFlagsVariable(varName)}
    elif KratosGlobals.Kernel.HasVariableData(varName):
        raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is defined but is of unsupported type\n".format(varName))
    else:
        raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is unknown. Maybe you need to import the application where it is defined?\n".format(varName))

def GetHistoricalVariableList(model_part, container, full_search=False):
    """ Get all variables in the historical database ( GetSolutionStepValues )

    Get all variables in the historical database ( GetSolutionStepValues ).
    If not specified the function will asume that all nodes in `model_part`
    have the same variables. 

    If `skip_full_check` is set to false a full
    check over all nodes will be done. This option can be slow.

    """

    variable = GetAsCommponents(GetVariableAndType(model_part.GetHistoricalVariables(container, full_search)))

    return []


def GetNonHistoricalVariableList(model_part):
    pass

def CheckAll(model_part):
    CheckAllHistoricalVariables(model_part)
    CheckAllNonHistoricalVariables(model_part)
    CheckAllVariablesStatistics(model_part)

def CheckAllHistoricalVariables(model_part):
    for variable in GetHistoricalVariableList(model_part):
        debug_utilities.CheckHistoricalVariable(model_part, variable)

def CheckAllNonHistoricalVariables(model_part):
    for variable in GetHistoricalVariableList(model_part, model_part.Nodes):
        debug_utilities.ChecknonHistoricalVariable(model_part, model_part.Nodes, variable)

    for variable in GetHistoricalVariableList(model_part, model_part.Elements):
        debug_utilities.ChecknonHistoricalVariable(model_part, model_part.Elements, variable)

    for variable in GetHistoricalVariableList(model_part, model_part.Conditions):
        debug_utilities.ChecknonHistoricalVariable(model_part, model_part.Conditions, variable)

def CheckAllVariablesStatistics(model_part):
    for variable in GetNonHistoricalVariableList(model_part):
        CheckVariableStatistics(model_part, variable)

def CheckVariableStatistics(model_part, input_variable):
    data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
    
    if data_comm.IsDistributed():
        mean_var = 0.0
        fixed_count = 0
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            fixed_count += int(node.IsFixed(input_variable))
            mean_var += node.GetSolutionStepValue(input_variable)

        global_mean_var = data_comm.SumAll(mean_var)
        global_nodes = data_comm.SumAll(model_part.GetCommunicator().LocalMesh().NumberOfNodes())
        global_fixed_count = data_comm.SumAll(fixed_count)

        print("Mean", input_variable.Name(), global_mean_var, global_nodes, global_mean_var/global_nodes, global_fixed_count)
    else:
        mean_var = 0.0
        fixed_count = 0
        for node in model_part.Nodes:
            fixed_count += int(node.IsFixed(input_variable))
            mean_var += node.GetSolutionStepValue(input_variable)

        print("Mean", input_variable.Name(), mean_var, model_part.NumberOfNodes(), mean_var/model_part.NumberOfNodes(), fixed_count) 