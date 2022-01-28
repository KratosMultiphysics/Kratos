# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(*args):
    return CreatePointLoadModelPart(*args)

class CreatePointLoadModelPart(CoSimulationCouplingOperation):
    """This operation creates a submodelpart containing PointLoad Conidtions for transferring loads
    """
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        IssueDeprecationWarning('CreatePointLoadModelPart', 'please use CreatePointBasedEntitiesProcess" instead')
        super().__init__(settings, process_info, data_communicator)
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
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "sub_model_part_name" : "UNSPECIFIED",
            "computing_model_part_name" : "UNSPECIFIED"
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults



