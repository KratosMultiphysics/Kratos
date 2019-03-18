from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *

def WriteFluidMdpa( filename, model_part, element_name, condition_name):
        restart_file = open(filename + ".mdpa", 'w')
        import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(model_part.Nodes, restart_file)
        new_restart_utilities.PrintElements(element_name,model_part.Elements, restart_file)
        new_restart_utilities.PrintConditions(condition_name,model_part.Conditions, restart_file)

        #printing info for slip BC
        new_restart_utilities.PrintRestart_Variable_On_Condition(IS_STRUCTURE, "IS_STRUCTURE", model_part.Conditions, restart_file)

        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_X,
            "VELOCITY_X",
            model_part.Nodes,
            restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_Y,
            "VELOCITY_Y",
            model_part.Nodes,
            restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_Z,
            "VELOCITY_Z",
            model_part.Nodes,
            restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            PRESSURE,
            "PRESSURE",
            model_part.Nodes,
            restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VISCOSITY,
            "VISCOSITY",
            model_part.Nodes,
            restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            DENSITY,
            "DENSITY",
            model_part.Nodes,
            restart_file)

        restart_file.close()
