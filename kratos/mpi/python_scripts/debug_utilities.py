import KratosMultiphysics
import KratosMultiphysics.mpi

def GetVariableAndType(var_name, split_by_commponents):
    """
    Return the variable with the given name.

    Keyword arguments:
    var_name -- The name of the variable to return
    split_by_commponents -- If set to True will split the variable into commponents if possible (VAR --> VAR_X, VAR_Y, etc...)
    """
    KratosGlobals = KratosMultiphysics.KratosGlobals

    if KratosGlobals.Kernel.HasDoubleVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetDoubleVariable(var_name)}
    elif KratosGlobals.Kernel.HasArrayVariable(var_name):
        if split_by_commponents:
            for component in ["X","Y","Z","XX","YY","XY"]:
                if KratosGlobals.Kernel.HasDoubleVariable(var_name+"_"+component):
                    yield {"name": var_name+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(var_name+"_"+component)}
        else:
            yield {"name": var_name, "value": KratosGlobals.Kernel.GetArrayVariable(var_name)}
    elif KratosGlobals.Kernel.HasArray4Variable(var_name):
        if split_by_commponents:
            for component in ["XX","XY","YX","YY"]:
                yield {"name": var_name+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(var_name+"_"+component)}
        else:
            yield {"name": var_name, "value": KratosGlobals.Kernel.GetArray4Variable(var_name)}
    elif KratosGlobals.Kernel.HasArray6Variable(var_name):
        if split_by_commponents:
            for component in ["XX","YY","ZZ","XY","YZ","XZ"]:
                yield {"name": var_name+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(var_name+"_"+component)}
        else:
            yield {"name": var_name, "value": KratosGlobals.Kernel.GetArray6Variable(var_name)}
    elif KratosGlobals.Kernel.HasArray9Variable(var_name):
        if split_by_commponents:
            for component in ["XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ"]:
                yield {"name": var_name+"_"+component, "value": KratosGlobals.Kernel.GetDoubleVariable(var_name+"_"+component)}
        else:
            yield {"name": var_name, "value": KratosGlobals.Kernel.GetArray9Variable(var_name)}
    elif KratosGlobals.Kernel.HasBoolVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetBoolVariable(var_name)}
    elif KratosGlobals.Kernel.HasIntVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetIntVariable(var_name)}
    elif KratosGlobals.Kernel.HasUnsignedIntVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetUnsignedIntVariable(var_name)}
    elif KratosGlobals.Kernel.HasVectorVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetVectorVariable(var_name)}
    elif KratosGlobals.Kernel.HasMatrixVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetMatrixVariable(var_name)}
    elif KratosGlobals.Kernel.HasStringVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetStringVariable(var_name)}
    elif KratosGlobals.Kernel.HasFlagsVariable(var_name):
        yield {"name": var_name, "value": KratosGlobals.Kernel.GetFlagsVariable(var_name)}
    elif KratosGlobals.Kernel.HasVariableData(var_name):
        raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is defined but is of unsupported type\n".format(var_name))
    else:
        raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is unknown. Maybe you need to import the application where it is defined?\n".format(var_name))

def GetHistoricalVariableList(model_part, split_by_commponents=False):
    """Get all variables in the historical database (GetSolutionStepValues)."""
    return GetVariableAndType(model_part.GetHistoricalVariableNames(), split_by_commponents)

def GetNodalNonHistoricalVariableList(modelPart, container, full_search=False, split_by_commponents=False):
    """
    Get all variables in the historical database (GetValues) for entities in `container`.

    Get all variables in the historical database (GetSolutionStepValues) for entities in `container`.
    If not specified the function will asume that all entities in `container` have the same variables.

    Keyword arguments:
    model_part -- Input Modelpart
    container -- Specific container of the modelpart to perform the check (Elements, Nodes, etc...)
    full_search -- Default: False. If set to true will individually check all the variables for every single individual entity in the container.
    split_by_commponents -- If set to True will split the variable into commponents if possible (VAR --> VAR_X, VAR_Y, etc...)
    """
    return GetVariableAndType(modelPart.GetNonHistoricalVariables(container, full_search), split_by_commponents)

def GetNonHistoricalVariableList(modelPart, full_search=False, split_by_commponents=False):
    """Get all variables in the non historical database (GetValues) for all entities(nodes, elements and conditions).

    Get all variables in the non historical database (GetValues) for all entities (nodes, elements and conditions).
    If not specified the function will asume that all entities of the same type in `modelPart` have the same variables.

    Keyword arguments:
    modelPart -- Input Modelpart
    full_search -- Default: False. If set to true will individually check all the variables for every single individual entity in the container.
    split_by_commponents -- If set to True will split the variable into commponents if possible (VAR --> VAR_X, VAR_Y, etc...)
    """
    node_variables = GetNodalNonHistoricalVariableList(modelPart.GetNonHistoricalVariables(modelPart.Nodes,      full_search, split_by_commponents))
    elem_variables = GetNodalNonHistoricalVariableList(modelPart.GetNonHistoricalVariables(modelPart.Elements,   full_search, split_by_commponents))
    cond_variables = GetNodalNonHistoricalVariableList(modelPart.GetNonHistoricalVariables(modelPart.Conditions, full_search, split_by_commponents))

    return node_variables + elem_variables + cond_variables

def CheckAll(model_part):
    """Check the consistency of all entities and its variables."""
    CheckAllHistoricalVariables(model_part)
    CheckAllNonHistoricalVariables(model_part)
    CheckAllVariablesStatistics(model_part)

def CheckAllHistoricalVariables(model_part):
    """Check the consistency of all historical variables."""
    for variable in GetHistoricalVariableList(model_part):
        debug_utilities.CheckHistoricalVariable(model_part, variable)

def CheckAllNonHistoricalVariables(model_part):
    """Check the consistency of all non-historical variables."""
    for variable in GetHistoricalVariableList(model_part, model_part.Nodes):
        debug_utilities.ChecknonHistoricalVariable(model_part, model_part.Nodes, variable)

    for variable in GetHistoricalVariableList(model_part, model_part.Elements):
        debug_utilities.ChecknonHistoricalVariable(model_part, model_part.Elements, variable)

    for variable in GetHistoricalVariableList(model_part, model_part.Conditions):
        debug_utilities.ChecknonHistoricalVariable(model_part, model_part.Conditions, variable)

def CheckAllHistoricalVariablesStatistics(model_part):
    """Compute the statistics of all the historical variables."""
    for variable in GetHistoricalVariableList(model_part):
        CheckVariableStatistics(model_part, variable)

def CheckAllNonHistoricalVariablesStatistics(model_part):
    """Compute the statistics of all the non-historical variables."""
    for variable in GetHistoricalVariableList(model_part, model_part.Nodes):
        CheckVariableStatistics(model_part, model_part.Nodes, variable)

    for variable in GetHistoricalVariableList(model_part, model_part.Elements):
        CheckVariableStatistics(model_part, model_part.Elements, variable)

    for variable in GetHistoricalVariableList(model_part, model_part.Conditions):
        CheckVariableStatistics(model_part, model_part.Conditions, variable)

def CheckHistoricalVariableStatistics(model_part, input_variable):
    """Compute the statistics of the given historical variable."""
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