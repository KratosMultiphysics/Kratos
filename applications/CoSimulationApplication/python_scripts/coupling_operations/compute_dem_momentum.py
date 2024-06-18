# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.DEMApplication import *

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation
import KratosMultiphysics.DEMFEMVolumeCouplingApplication as VCA

from collections import deque


def Create(*args):
    return ComputeDemMomentum(*args)


class ComputeDemMomentum(CoSimulationCouplingOperation):

    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = self.model[self.model_part_name]
        self.force_end_time=0.025
        self.dt =  1e-5
        self.step=0
        self.timesteps=self.force_end_time/self.dt
        self.max_force=0.025
        self.force_slope=self.max_force/self.timesteps # to reach max force within force_end_time
        #self.displacements_history = {node.Id: deque(maxlen=1000) for node in self.model_part.Nodes}

    def InitializeCouplingIteration(self):
        
        utils = VCA.DEMFEMVolumeCouplingUtilities()

        utils.CalculateMomentum(self.model_part)
        # for node in self.model_part.Nodes:
        #     #current_displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
            
        #     # Append current displacement to the deque
        #     #self.displacements_history[node.Id].append(current_displacement)
            
        #     # Get displacement from 1000 timesteps ago, if available
        #     #displacement_1000_timesteps_ago = self.displacements_history[node.Id][0] if len(self.displacements_history[node.Id]) == 1000 else None
            
        #     # Calculate the difference in displacement
        #     # if displacement_1000_timesteps_ago:
        #     #     displacement_new= current_displacement - displacement_1000_timesteps_ago
                
        #     #     # Use this displacement difference for your calculations
        #     #     node.SetSolutionStepValue(VCA.DISPLACEMENT_MULTIPLIED_MASS, node.GetSolutionStepValue(KM.NODAL_MASS) * displacement_new)
        #     # else:
        #     #     # If there's no displacement from 1000 timesteps ago (like in the first few timesteps), you can choose to do nothing or use the current displacement
        #     node.SetSolutionStepValue(VCA.DISPLACEMENT_MULTIPLIED_MASS, node.GetSolutionStepValue(KM.NODAL_MASS)* node.GetSolutionStepValue(KM.VELOCITY))

    def FinalizeCouplingIteration(self):

        utils = VCA.DEMFEMVolumeCouplingUtilities()

        utils.CalculateDEMForces(self.model_part)
        # node_ids = [22,40,47,52] # for assigning loads to the bottom nodes
        # self.step+=1
        # self.force_end_time=0.025
        # self.dt=  1e-5
        # self.timesteps=self.force_end_time/self.dt
        # self.max_force=0.025
        # #self.max_force=0.00625 # for bigger domain
        # self.force_slope=self.max_force/self.timesteps
        # pointload=min(self.force_slope*self.step,self.max_force)
        for node in self.model_part.Nodes:

            print("#############################################      node.Id ",node.Id," COUPLING FORCE ON PARTICLE  = ",node.GetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE))
            
        #     particle_weight = node.GetSolutionStepValue(VCA.PARTICLE_COUPLING_WEIGHT)
        #     if(particle_weight==0):
        #         node.SetSolutionStepValue(VCA.PARTICLE_COUPLING_WEIGHT,1)
        #     else:
        #         node.SetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE, node.GetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE)*node.GetSolutionStepValue(KM.NODAL_MASS)) # if nodal mass needs to be multiplied with external applied force
            # print("#############################################      node.Id ",node.Id," PARTICLE_COUPLING_WEIGHT = ",node.GetSolutionStepValue(VCA.PARTICLE_COUPLING_WEIGHT))
            # if node.Id in node_ids:
            #     node.SetSolutionStepValue(KM.EXTERNAL_APPLIED_FORCE, [0,pointload,0])
    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : ""
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
