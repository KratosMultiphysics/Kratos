from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

import assign_vector_components_to_nodes_process as BaseProcess

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignInletProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignInletProcess(BaseProcess.AssignVectorComponentsToNodesProcess):
    def __init__(self, Model, custom_settings ):

        super(AssignInletProcess, self).__init__(Model, custom_settings)

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        entity_type = "Nodes"
        assign_flags = [KratosMultiphysics.INLET]
        KratosSolid.AssignFlagsToEntitiesProcess(self.model_part,entity_type,assign_flags).Execute()

        super(AssignInletProcess, self).ExecuteInitialize()

        for process in self.AssignValueProcesses:
            process.ExecuteInitialAssignment()
