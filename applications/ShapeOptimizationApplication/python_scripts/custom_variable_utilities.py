# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================
def WriteDictionaryDataOnNodalVariable(data, mode_part, nodal_variable):
    for node_id, tmp_gradient in data.items():
        mode_part.Nodes[node_id].SetSolutionStepValue(nodal_variable,0,tmp_gradient)

# ------------------------------------------------------------------------------
def ReadNodalVariableToList(mode_part, nodal_variable):
    variable_values_list = []
    for node in mode_part.Nodes:
        tmp_values = node.GetSolutionStepValue(nodal_variable)
        variable_values_list.append(tmp_values[0])
        variable_values_list.append(tmp_values[1])
        variable_values_list.append(tmp_values[2])
    return variable_values_list

# ==============================================================================
