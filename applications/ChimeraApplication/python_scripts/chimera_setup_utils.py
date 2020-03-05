from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera


def SeparateAndValidateChimeraSettings(parameters):
    # Deprecation warnings
    solver_parameters = parameters
    # Checking if the parameters has 'chimera_settings' entry.
    # This is required.
    if solver_parameters.Has("chimera_settings"):
        chimera_parameters = solver_parameters["chimera_settings"].Clone()
    else:
        raise Exception("The \"solver_settings\" should have the entry \"chimera_settings\" ")

    # Seperating the fluid solver settings.
    if solver_parameters.Has("fluid_solver_settings"):
        fluid_parameters = solver_parameters["fluid_solver_settings"].Clone()
    else:
        fluid_parameters = solver_parameters.Clone()

    # Extracting the chimera_parts. this is required for ApplyChimera process.
    if chimera_parameters.Has("chimera_parts"):
        chimera_levels = chimera_parameters["chimera_parts"].Clone()
    else:
        raise Exception("The \"solver_settings\" should have the entry \"chimera_parts\" ")

    chimera_echo_lvl = 0
    if chimera_parameters.Has("chimera_echo_level"):
        chimera_echo_lvl = chimera_parameters["chimera_echo_level"].GetInt()
    else:
        chimera_echo_lvl = fluid_parameters["echo_level"].GetInt()

    chimera_internal_parts = []
    for level in chimera_levels:
        for level_parameters in level :
            if level_parameters.Has("internal_parts_for_chimera"):
                part_name_list = level_parameters["internal_parts_for_chimera"]
                for part_name in part_name_list:
                    chimera_internal_parts.append(part_name.GetString())

    reformulate_every_step = False
    if chimera_parameters.Has("reformulate_chimera_every_step"):
        reformulate_every_step = chimera_parameters["reformulate_chimera_every_step"].GetBool()

    if not fluid_parameters.Has("reform_dofs_at_each_step"):
        fluid_parameters.AddEmptyValue("reform_dofs_at_each_step")

    fluid_parameters["reform_dofs_at_each_step"].SetBool(reformulate_every_step)
    solver_parameters.RemoveValue("chimera_settings")

    return [chimera_parameters, chimera_internal_parts, parameters]

def GetApplyChimeraProcess(model, chimera_parameters, fluid_parameters):

    chimera_echo_lvl = 0
    if chimera_parameters.Has("chimera_echo_level"):
        chimera_echo_lvl = chimera_parameters["chimera_echo_level"].GetInt()
    else:
        chimera_echo_lvl = fluid_parameters["echo_level"].GetInt()

    reformulate_every_step = False
    if chimera_parameters.Has("reformulate_chimera_every_step"):
        reformulate_every_step = chimera_parameters["reformulate_chimera_every_step"].GetBool()

    # Extracting the chimera_parts. this is required for ApplyChimera process.
    if chimera_parameters.Has("chimera_parts"):
        chimera_levels = chimera_parameters["chimera_parts"].Clone()
    else:
        raise Exception("The \"solver_settings\" should have the entry \"chimera_parts\" ")

    main_model_part = model[fluid_parameters["model_part_name"].GetString()]
    domain_size = main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
    solver_type = fluid_parameters["solver_type"].GetString()


    # Creating the necessary variant of the apply chimera process.
    if domain_size == 2:
        if(solver_type == "Monolithic" or solver_type == "monolithic"):
            chimera_process = KratosChimera.ApplyChimeraProcessMonolithic2d(main_model_part,chimera_levels)
        elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
            chimera_process = KratosChimera.ApplyChimeraProcessFractionalStep2d(main_model_part,chimera_levels)
    else:
        if(solver_type == "Monolithic" or solver_type == "monolithic"):
            chimera_process = KratosChimera.ApplyChimeraProcessMonolithic3d(main_model_part,chimera_levels)
        elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
            chimera_process = KratosChimera.ApplyChimeraProcessFractionalStep3d(main_model_part,chimera_levels)

    chimera_process.SetEchoLevel(chimera_echo_lvl)
    chimera_process.SetReformulateEveryStep(reformulate_every_step)

    return chimera_process


def SetChimeraInternalPartsFlag(model, chimera_internal_parts):
    '''
        This function sects the flag CHIMERA_INTERNAL_BOUNDARY on the specified modelparts
        so that they are excluded from the extract surface operation later on.
    '''
    for mp_name in chimera_internal_parts:
        KratosMultiphysics.VariableUtils().SetFlag(KratosChimera.CHIMERA_INTERNAL_BOUNDARY, True,  model[mp_name].Nodes)