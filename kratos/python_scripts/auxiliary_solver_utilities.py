from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
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
                dof_variable_x = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name + "_X")
                reaction_variable_x = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name + "_X")
                KratosMultiphysics.VariableUtils().AddDof(dof_variable_x, reaction_variable_x, model_part)
                dof_variable_y = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name + "_Y")
                reaction_variable_y = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name + "_Y")
                KratosMultiphysics.VariableUtils().AddDof(dof_variable_y, reaction_variable_y, model_part)
                dof_variable_z = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name + "_Z")
                reaction_variable_z = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name + "_Z")
                KratosMultiphysics.VariableUtils().AddDof(dof_variable_z, reaction_variable_z, model_part)
            else:
                raise Exception("The variable " + dof_variable_name + " is incompatible in type")
        else:
            err_mssg = "The added reaction variable: " + reaction_variable_name
            err_mssg+= "(" + KratosMultiphysics.KratosGlobals.GetVariableType(reaction_variable_name) + ") "
            err_mssg+= "is not the same type as the pair DoF: " + dof_variable_name
            err_mssg+= "(" + KratosMultiphysics.KratosGlobals.GetVariableType(dof_variable_name) + ")"
            raise Exception(err_mssg)


