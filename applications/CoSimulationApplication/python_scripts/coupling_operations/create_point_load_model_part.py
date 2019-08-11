from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrappers):
    cs_tools.SettingsTypeCheck(settings)
    return CreatePointLoadModelPart(settings, solver_wrappers)

class CreatePointLoadModelPart(CoSimulationCouplingOperation):
    """This operation creates a submodelpart containing PointLoad Conidtions for transferring loads
    """
    def __init__(self, settings, solver_wrappers):
        super(CreatePointLoadModelPart, self).__init__(settings)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model

    def Initialize(self):
        computing_model_part_name = self.settings["computing_model_part_name"].GetString()
        sub_model_part_name = self.settings["sub_model_part_name"].GetString()

        computing_domain = self.model[computing_model_part_name]

        node_id_list = []
        number_of_conditions = computing_domain.NumberOfConditions()
        for cond_counter,node_i in enumerate(computing_domain.Nodes):
            node_id_list.append(node_i.Id)
            computing_domain.CreateNewCondition(
                "PointLoadCondition3D1N",cond_counter+number_of_conditions+1,[node_i.Id],
                computing_domain.GetProperties()[0])

        struct_smp = computing_domain.CreateSubModelPart(sub_model_part_name)
        struct_smp.AddNodes(node_id_list)

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "sub_model_part_name" : "UNSPECIFIED",
            "computing_model_part_name" : "UNSPECIFIED"
        }""")
        this_defaults.AddMissingParameters(super(CreatePointLoadModelPart, cls)._GetDefaultSettings())
        return this_defaults



