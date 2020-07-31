import KratosMultiphysics
import KratosMultiphysics.mpi
import KratosMultiphysics.mpi.debug_utilities as debug_utilities

def GetVariable(self, VarName):
    """ This method returns the variable with the given name

    Keyword arguments:
    self -- It signifies an instance of a class.
    VarName -- The name of the variable to return

    """

    if KratosGlobals.Kernel.HasDoubleVariable(VarName):
        return KratosGlobals.Kernel.GetDoubleVariable(VarName)
    elif KratosGlobals.Kernel.HasArrayVariable(VarName):
        return KratosGlobals.Kernel.GetArrayVariable(VarName)
    elif KratosGlobals.Kernel.HasArray4Variable(VarName):
        return KratosGlobals.Kernel.GetArray4Variable(VarName)
    elif KratosGlobals.Kernel.HasArray6Variable(VarName):
        return KratosGlobals.Kernel.GetArray6Variable(VarName)
    elif KratosGlobals.Kernel.HasArray9Variable(VarName):
        return KratosGlobals.Kernel.GetArray9Variable(VarName)
    elif KratosGlobals.Kernel.HasBoolVariable(VarName):
        return KratosGlobals.Kernel.GetBoolVariable(VarName)
    elif KratosGlobals.Kernel.HasIntVariable(VarName):
        return KratosGlobals.Kernel.GetIntVariable(VarName)
    elif KratosGlobals.Kernel.HasUnsignedIntVariable(VarName):
        return KratosGlobals.Kernel.GetUnsignedIntVariable(VarName)
    elif KratosGlobals.Kernel.HasVectorVariable(VarName):
        return KratosGlobals.Kernel.GetVectorVariable(VarName)
    elif KratosGlobals.Kernel.HasMatrixVariable(VarName):
        return KratosGlobals.Kernel.GetMatrixVariable(VarName)
    elif KratosGlobals.Kernel.HasStringVariable(VarName):
        return KratosGlobals.Kernel.GetStringVariable(VarName)
    elif KratosGlobals.Kernel.HasFlagsVariable(VarName):
        return KratosGlobals.Kernel.GetFlagsVariable(VarName)
    elif KratosGlobals.Kernel.HasVariableData(VarName):
        raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is defined but is of unsupported type\n".format(VarName))
    else:
        raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is unknown. Maybe you need to import the application where it is defined?\n".format(VarName))

def GetHistoricalVariableList(model_part, skip_full_check=True):
    """ Get all variables in the historical database ( GetSolutionStepValues )

    Get all variables in the historical database ( GetSolutionStepValues ).
    If not specified the function will asume that all nodes in `model_part`
    have the same variables. 
    
    If `skip_full_check` is set to false a full
    check over all nodes will be done. This option can be slow.

    """

    variable_names = KratosGlobals.Kernel.GetallVariablesNames()

    # if skip_full_check:
    # else:
    #     for node in model_part.Nodes:
    #         for name in variable_names:
    #             yield GetVariable(name)


def GetNonHistoricalVariableList(model_part):

def CheckAll(model_part):
    CheckAllHistoricalVariables(model_part)
    CheckAllNonHistoricalVariables(model_part)
    CheckAllVariablesStatistics(model_part)

def CheckAllHistoricalVariables(model_part):
    for variable in GetHistoricalVariableList(model_part):
        debug_utilities.CheckHistoricalVariable(model_part, variable)

def CheckAllNonHistoricalVariables(model_part):
    for variable in GetNonHistoricalVariableList(model_part):
        debug_utilities.CheckNonHistoricalVariable(model_part, variable)

def CheckAllVariablesStatistics():
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

        print("Mean", input_variable.Name(), mean_var, model_part.NumberOfNodes(), mean_var/model_part.NumberOfNodes(), fixed_count)% 