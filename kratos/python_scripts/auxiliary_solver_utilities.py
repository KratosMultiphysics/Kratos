import KratosMultiphysics

def AddVariables(model_part, aux_variable_list):
  for variable_name in aux_variable_list.GetStringArray():
      variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
      model_part.AddNodalSolutionStepVariable(variable)

def AddDofs(model_part, aux_dofs_list, aux_reaction_list):
    if (aux_dofs_list.size() != aux_reaction_list.size()):
        raise Exception("DoFs list and reaction list should be the same size")
    for dof_variable_name, reaction_variable_name in  zip(aux_dofs_list.GetStringArray(), aux_reaction_list.GetStringArray()):
        dof_variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(dof_variable_name)
        reaction_variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(reaction_variable_name)
        if(dof_variable_type == reaction_variable_type):
            if(dof_variable_type == "Double"): # Double variable
                dof_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name)
                reaction_variable = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name)
                KratosMultiphysics.VariableUtils().AddDof(dof_variable, reaction_variable, model_part)
            elif (dof_variable_type == "Array"): # Components variable
                for comp in ["_X", "_Y", "_Z"]:
                    dof_variable_comp = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name + comp)
                    reaction_variable_comp = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name + comp)
                    KratosMultiphysics.VariableUtils().AddDof(dof_variable_comp, reaction_variable_comp, model_part)
            else:
                raise Exception("The variable " + dof_variable_name + " is incompatible in type")
        else:
            err_msg = "The added reaction variable: " + reaction_variable_name
            err_msg+= "(" + KratosMultiphysics.KratosGlobals.GetVariableType(reaction_variable_name) + ") "
            err_msg+= "is not the same type as the pair DoF: " + dof_variable_name
            err_msg+= "(" + KratosMultiphysics.KratosGlobals.GetVariableType(dof_variable_name) + ")"
            raise Exception(err_msg)


