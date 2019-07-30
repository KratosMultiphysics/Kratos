from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics

def add_auxiliary_variables(model_part, aux_variable_list):
  for i in range(aux_variable_list.size()):
      variable_name = aux_variable_list[i].GetString()
      variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
      model_part.AddNodalSolutionStepVariable(variable)

def add_auxiliary_dofs(model_part, aux_dofs_list, aux_reaction_list):
  if (aux_dofs_list.size() != aux_reaction_list.size()):
          raise Exception("DoFs list and reaction list should be the same size")
  for i in range(aux_dofs_list.size()):
      dof_variable_name = aux_dofs_list[i].GetString()
      reaction_variable_name = aux_reaction_list[i].GetString()
      if (KratosMultiphysics.KratosGlobals.HasVariable(dof_variable_name)):
          if(KratosMultiphysics.KratosGlobals.GetVariableType(dof_variable_name) == "Double"): # Double variable
              dof_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name)
              reaction_variable = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name)
              KratosMultiphysics.VariableUtils().AddDof(dof_variable, reaction_variable, model_part)
          elif (KratosMultiphysics.KratosGlobals.GetVariableType(dof_variable_name) == "Array"): # Components variable
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
          KratosMultiphysics.Logger.PrintWarning("auxiliary_reaction_list list", "The variable " + dof_variable_name + "is not a compatible type")


