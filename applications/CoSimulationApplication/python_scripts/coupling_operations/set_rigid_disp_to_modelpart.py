from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ChimeraApplication as ChimeraApp

import numpy as np
import math

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrappers):
    cs_tools.SettingsTypeCheck(settings)
    return SetRigidDisplacementOperation(settings, solver_wrappers)

class SetRigidDisplacementOperation(CoSimulationCouplingOperation):
    """This operation sets the Mesh displacement of a given master node
        to the specified modelpart. This will also do the following operation
        on the nodes of the modelpart : x = x0+disp_of_master_node
        and sets velocity of master node on to the modelpart
    """
    def __init__(self, settings, solver_wrappers):
        super(SetRigidDisplacementOperation, self).__init__(settings)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeCouplingIteration(self):
        pass

    def FinalizeCouplingIteration(self):
        pass

    def Execute(self):
        self.modelpart = self.interface_data.GetModelPart()
        master_node_id = self.settings["master_node_id"].GetInteger()
        master_node = self.modelpart.Nodes[master_node_id]
        master_node_x_disp = master_node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X)
        master_node_y_disp = master_node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y)
        master_node_z_disp = master_node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Z)
        for node in self.modelpart.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X, 0, master_node_x_disp)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y, 0, master_node_y_disp)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Z, 0, master_node_z_disp)

            node.X = node.X0 + master_node_x_disp
            node.Y = node.Y0 + master_node_y_disp
            node.Z = node.Z0 + master_node_z_disp

            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, 0, master_node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY))

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "master_node_id" : 0
        }""")
        this_defaults.AddMissingParameters(super(SetRigidDisplacementOperation, cls)._GetDefaultSettings())
        return this_defaults
