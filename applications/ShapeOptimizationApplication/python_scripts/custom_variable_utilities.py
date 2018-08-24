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
    data_dimension = len(next(iter(data.values())))
    if data_dimension != 3:
        raise RuntimeError("custom_variable_utility::WriteDictionaryDataOnNodalVariable: Wrong dimension of data entries in input dictionary!")

    for node_id, tmp_gradient in data.items():
        model_part.Nodes[node_id].SetSolutionStepValue(nodal_variable, 0, tmp_gradient)

# ------------------------------------------------------------------------------
def ReadNodalVariableToList(model_part, nodal_variable):
    variable_values_list = [0.0]*model_part.NumberOfNodes()*3

    for itr, node in enumerate(model_part.Nodes):
        tmp_values = node.GetSolutionStepValue(nodal_variable)
        variable_values_list[3*itr+0] = tmp_values[0]
        variable_values_list[3*itr+1] = tmp_values[1]
        variable_values_list[3*itr+2] = tmp_values[2]

    return variable_values_list

# --------------------------------------------------------------------------
def WriteListToNodalVariable(input_list, model_part, nodal_variable):
    if len(input_list) != model_part.NumberOfNodes()*3:
        raise RuntimeError("custom_variable_utility::WriteListToNodalVariable: Wrong size of input_list!")

    k = 0
    for node in model_part.Nodes:
        tmp_values = input_list[k:k+3]
        node.SetSolutionStepValue(nodal_variable, tmp_values)
        k = k+3

# --------------------------------------------------------------------------
def WriteNodeCoordinatesToList(model_part):
    variable_values_list = [0.0]*model_part.NumberOfNodes()*3

    for itr, node in enumerate(model_part.Nodes):
        variable_values_list[3*itr+0] = node.X
        variable_values_list[3*itr+1] = node.Y
        variable_values_list[3*itr+2] = node.Z

    return variable_values_list

# ==============================================================================
