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

def AddAuxiliaryDofsToDofsWithReactionsList(aux_dofs_list, aux_reaction_list, dofs_with_reactions_list):
    if not aux_dofs_list.size() == aux_reaction_list.size():
        raise Exception("DOFs list and reaction list should be the same size")
    for dof_variable_name, reaction_variable_name in  zip(aux_dofs_list.GetStringArray(), aux_reaction_list.GetStringArray()):
        dof_variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(dof_variable_name)
        reaction_variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(reaction_variable_name)
        if dof_variable_type == reaction_variable_type:
            if dof_variable_type == "Double":
                dofs_with_reactions_list.append([dof_variable_name, reaction_variable_name])
            elif dof_variable_type == "Array":
                for comp in ["_X", "_Y", "_Z"]:
                    dofs_with_reactions_list.append([dof_variable_name + comp, reaction_variable_name + comp])
            else:
                raise Exception("The variable " + dof_variable_name + " is incompatible in type")
        else:
            err_msg = "The added reaction variable \'{}\' (\'{}\') is not the same type as the pair DOF \'{}\' (\'{}\').".format(reaction_variable_name, reaction_variable_type, dof_variable_name, dof_variable_type)
            raise Exception(err_msg)

def SetAndFillBuffer(model_part, required_buffer_size, delta_time):
    # Check input data
    if required_buffer_size < 1:
        raise Exception("Provided buffer size is {}. Value greater or equal to one must be provided.".format(required_buffer_size))
    if delta_time < 0.0:
        raise Exception("Provided delta time is {}. A value greater or equal to zero must be provided.".format(delta_time))

    # Set the buffer size for the nodal solution steps data.
    # Existing nodal solution step data may be lost.
    current_buffer_size = model_part.GetBufferSize()
    buffer_size = max(current_buffer_size, required_buffer_size)
    model_part.SetBufferSize(buffer_size)

    # Cycle the buffer. This sets all historical nodal solution step data to
    # the current value and initializes the time stepping in the process info.
    step = -buffer_size
    init_time = model_part.ProcessInfo[KratosMultiphysics.TIME]
    time = init_time - delta_time * buffer_size
    model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
    for _ in range(buffer_size):
        step += 1
        time += delta_time
        model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
        model_part.CloneTimeStep(time)
