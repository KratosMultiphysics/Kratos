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
def ReadNodalVariableToList(model_part, nodal_variable):
    variable_values_list = []
    for node in model_part.Nodes:
        tmp_values = node.GetSolutionStepValue(nodal_variable)
        variable_values_list.append(tmp_values[0])
        variable_values_list.append(tmp_values[1])
        variable_values_list.append(tmp_values[2])
    return variable_values_list

# --------------------------------------------------------------------------
def WriteListToNodalVariable(input_list, model_part, nodal_variable):
    k = 0
    for node in model_part.Nodes:
        tmp_values = input_list[k:k+3]
        node.SetSolutionStepValue(nodal_variable, tmp_values)
        k = k+3

# --------------------------------------------------------------------------
def WriteNodeCoordinatesToList(model_part):
    list_of_values = []
    for node in model_part.Nodes:
        list_of_values.append(node.X)
        list_of_values.append(node.Y)
        list_of_values.append(node.Z)
    return list_of_values

# ==============================================================================
