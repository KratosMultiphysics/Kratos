# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================
def WriteDictionaryDataOnNodalVariable(data, model_part, nodal_variable):
    for node_id, tmp_gradient in data.items():
        model_part.Nodes[node_id].SetSolutionStepValue(nodal_variable, 0, tmp_gradient)

# ------------------------------------------------------------------------------
def ReadNodalVariableToList(model_part, nodal_variable, dimension=3):
    variable_values_list = [0.0]*model_part.NumberOfNodes()*dimension

    for itr, node in enumerate(model_part.Nodes):
        tmp_values = node.GetSolutionStepValue(nodal_variable)

        if dimension == 1:
            variable_values_list[itr] = tmp_values
        else:
            for i in range(dimension):
                variable_values_list[dimension*itr+i] = tmp_values[i]

    return variable_values_list

# --------------------------------------------------------------------------
def WriteListToNodalVariable(input_list, model_part, nodal_variable, dimension=3):
    if len(input_list) != model_part.NumberOfNodes()*dimension:
        raise RuntimeError("custom_variable_utility::WriteListToNodalVariable: Wrong size of input_list!")

    k = 0
    for node in model_part.Nodes:
        if dimension == 1:
            tmp_values = input_list[k]
        else:
            tmp_values = input_list[k:k+dimension]
        node.SetSolutionStepValue(nodal_variable, tmp_values)
        k = k+dimension

# --------------------------------------------------------------------------
def WriteNodeCoordinatesToList(model_part):
    variable_values_list = [0.0]*model_part.NumberOfNodes()*3

    for itr, node in enumerate(model_part.Nodes):
        variable_values_list[3*itr+0] = node.X
        variable_values_list[3*itr+1] = node.Y
        variable_values_list[3*itr+2] = node.Z

    return variable_values_list

# ==============================================================================
