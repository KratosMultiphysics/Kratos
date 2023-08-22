# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.DEMApplication import *

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation
import KratosMultiphysics.DEMFEMVolumeCouplingApplication as VCA

def Create(*args):
    return ComputeDemMomentum(*args)

class ComputeDemMomentum(CoSimulationCouplingOperation):

    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = self.model[self.model_part_name]
  

    def InitializeCouplingIteration(self): # currently total lagrange method , linear momentum contains mass x displacement
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(VCA.DISPLACEMENT_MULTIPLIED_MASS, node.GetSolutionStepValue(KM.NODAL_MASS)* node.GetSolutionStepValue(KM.DISPLACEMENT))
            #print("For node id:",node.Id,", PARTICLE_COUPLING_WEIGHT=",node.GetSolutionStepValue(VCA.PARTICLE_COUPLING_WEIGHT))
            #print("For node id:",node.Id,", PARTICLE_COUPLING_FORCE=",node.GetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE))
            #print("For node id:",node.Id,", EXTERNAL_APPLIED_FORCE=",node.GetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE))
            
    def FinalizeCouplingIteration(self):
        #node_ids = [22,40,47,52] # for assigning loads to the bottom nodes
        for node in self.model_part.Nodes:
            particle_weight = node.GetSolutionStepValue(VCA.PARTICLE_COUPLING_WEIGHT)
            if particle_weight!=0:
                node.SetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE, node.GetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE)/ particle_weight) # if nodal mass needs to be multiplied with external applied force
            #if node.Id in node_ids:
                #node.SetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE, [0,0.025,0])
    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : ""
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
