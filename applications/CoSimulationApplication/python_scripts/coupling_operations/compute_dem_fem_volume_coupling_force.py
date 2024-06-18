# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA 
import KratosMultiphysics.DEMFEMVolumeCouplingApplication as VCA

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

import math
from collections import deque
def Create(*args):
    return ComputeNodalCouplingForce(*args)

class ComputeNodalCouplingForce(CoSimulationCouplingOperation):

    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = self.model[self.model_part_name]
        self.penalty_max = self.settings["penalty_max"].GetDouble()  
        self.dt= self.settings["timestep"].GetDouble() 
        self.force_end_time= self.settings["force_end_time"].GetDouble() 
        self.tolerance= self.settings["boundary_tolerance"].GetDouble() 
        self.y_fem_boundary= self.settings["y_fem_boundary"].GetDouble() 
        self.y_dem_boundary= self.settings["y_dem_boundary"].GetDouble() 
        self.weight_fem_boundary = self.settings["weight_fem_boundary"].GetDouble() 
        self.weight_dem_boundary = self.settings["weight_dem_boundary"].GetDouble() 
        self.step=0
        self.timesteps=self.force_end_time/self.dt
        #self.max_force=0.025
        self.max_force=0.025
        #self.max_force=0.00625 # for bigger domain
        self.force_slope=self.max_force/self.timesteps # to reach max force within force_end_time
        #self.displacements_history = {node.Id: deque(maxlen=1000) for node in self.model_part.Nodes}


    def InitializeCouplingIteration(self):


        self.step+=1
        
        pointload=-1*min(self.force_slope*self.step,self.max_force) # for 4 nodes gradually increasing point load
        #pointload=-1*4*(1/9)*min(self.force_slope*self.step,self.max_force) # for 9 nodes gradually increasing point load

        node_ids = [1,2,3,4]  # for assigning point loads to the top 4 nodes
        #node_ids =[145,146,147,148,149,150,151,152,153]  # for assigning point loads to the top 9 nodes
        #node_ids = [1,2,3,7,10,11,15,17,23,27,28,30,33,41,43,60]  # for assigning point loads to the top nodes


        utils = VCA.DEMFEMVolumeCouplingUtilities()

        utils.AssignPointLoads(self.model_part,node_ids,pointload) # assigning point loads to the top nodes
        utils.SetNodalCouplingWeightsOnFEMLinearly(self.model_part,self.y_fem_boundary,self.y_dem_boundary,self.tolerance,self.weight_fem_boundary,self.weight_dem_boundary) 
        utils.CalculateDisplacementDifference(self.model_part,self.dt) # calculating displacement difference
        utils.CalculateNodalCouplingForces(self.model_part,self.penalty_max) # calculating nodal coupling forces including point load
        utils.CalculateNodalDEMCouplingForces(self.model_part)


        # for node in self.model_part.Nodes:  
          
        #         if node.Id in node_ids:   # assigning point loads - only if force is exerted on FEM.

        #             node.SetSolutionStepValue(SMA.POINT_LOAD,[0,pointload,0]) # for 4 nodes
                    #print("For node id:",node.Id,", point load=",node.GetSolutionStepValue(SMA.POINT_LOAD))

                #current_displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
            
                # Append current displacement to the deque
                #self.displacements_history[node.Id].append(current_displacement)
                
                # Get displacement from 1000 timesteps ago, if available
                #displacement_1000_timesteps_ago = self.displacements_history[node.Id][0] if len(self.displacements_history[node.Id]) == 1000 else None
                
                
                # tolerance = 0.001 # tolerance on the boundary limits
                # if node.Y >= self.y_fem_boundary - tolerance and node.Y <= self.y_dem_boundary + tolerance: # assigning weights to the nodes in the hybrid region
                #     weight= self.weight_fem_boundary + (self.weight_dem_boundary-self.weight_fem_boundary)*(node.Y-self.y_fem_boundary)/(self.y_dem_boundary-self.y_fem_boundary)
                #     #weight=0.5 #putting 0.5 as weight for testing
                #     node.SetSolutionStepValue(VCA.NODAL_COUPLING_WEIGHT, weight)
                #     #print("iffffffffffffffff")
                # elif node.Y <= self.y_fem_boundary + tolerance and node.Y >= self.y_dem_boundary - tolerance: # assigning weights to the nodes in the hybrid region
                #     weight= self.weight_dem_boundary + (self.weight_fem_boundary-self.weight_dem_boundary)*(node.Y-self.y_dem_boundary)/(self.y_fem_boundary-self.y_dem_boundary)
                #     #weight=0.5 #putting 0.5 as weight for testing
                #     node.SetSolutionStepValue(VCA.NODAL_COUPLING_WEIGHT, weight)
                #     #print("eeeellllllliiiiiiiffffffffffff")
                # else: 
                #     node.SetSolutionStepValue(VCA.NODAL_COUPLING_WEIGHT, 0) # assigning weights to the nodes in the non-hybrid region
                #     #print("elseeeeeeeeeeeeeeeeeeeeeeee")
               
                # total_mass = node.GetSolutionStepValue(KM.NODAL_MAUX)

                # if(total_mass!=0): # check for hybrid region
                #     node.SetSolutionStepValue(SMA.POINT_LOAD,[0,0,0]) #nodal coupling forces set to zero before every iteration
                #     node.SetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE,[0,0,0]) # force to map to dem (set to zero before every iteration)
                #     Velocity_dem = node.GetSolutionStepValue(VCA.DISPLACEMENT_MULTIPLIED_MASS) / total_mass # updated lagrange method-> calculating homogenised velocity
                #     # if displacement_1000_timesteps_ago:
                #     #     displacement_new = current_displacement - displacement_1000_timesteps_ago
                #     #     node.SetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT, Displacement_dem - displacement_new)# total lagrange method : calcualting displacement difference
                #     # else:
                #     node.SetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT,node.GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT)+ (Velocity_dem - node.GetSolutionStepValue(KM.VELOCITY))*self.dt)# updated lagrange method : calcualting displacement difference

        # for elem in self.model_part.Elements:
        #      if(elem.GetNodes()[0].GetSolutionStepValue(KM.NODAL_MAUX))!=0: 
        #         #V=[0,0,0,0,0,0,0,0]
        #         V=elem.CalculateOnIntegrationPoints(KM.INTEGRATION_WEIGHT,self.model_part.ProcessInfo) 
        #         for i in range(elem.GetGeometry().IntegrationPointsNumber()): #gauss quadrature 
        #             J = (elem.GetGeometry().DeterminantOfJacobian(i))
        #             #get integration weight at integration points
        #             #w = V[i]/J
        #             w = 1
        #             #print('integration weight', w/J)
        #             shape_functions = elem.GetGeometry().ShapeFunctionsValues()
        #             for n in range(len(elem.GetNodes())):
        #                 for m in range(len(elem.GetNodes())):
        #                     vol = self.penalty_max * w * J * shape_functions[(i,n)] * shape_functions[(i,m)]
        #                     elem.GetNodes()[n].SetSolutionStepValue(SMA.POINT_LOAD, elem.GetNodes()[n].GetSolutionStepValue(SMA.POINT_LOAD) + 1* vol * elem.GetNodes()[n].GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT))
        #                     elem.GetNodes()[n].SetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE,elem.GetNodes()[n].GetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE) -1* vol * elem.GetNodes()[n].GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT)) #storing force to map to dem 
        #print("After calculation of point loads")
        # for node in self.model_part.Nodes: 
        #     if(node.GetSolutionStepValue(KM.NODAL_MAUX))!=0:
                #print("For node id:",node.Id,", for particle coupling force=",node.GetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE))
                #node.SetSolutionStepValue(SMA.POINT_LOAD, node.GetSolutionStepValue(SMA.POINT_LOAD) / (1-node.GetSolutionStepValue(VCA.NODAL_COUPLING_WEIGHT))) # dividing by nodal coupling weight, replace 0 by VCA.NODAL_COUPLING_WEIGHT in the future.
                #  print("For node id:",node.Id,", point load=",node.GetSolutionStepValue(SMA.POINT_LOAD))
                #print("For node id:",node.Id,", nodal coupling weight=",node.GetSolutionStepValue(VCA.NODAL_COUPLING_WEIGHT))
                #node.SetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE, node.GetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE) / node.GetSolutionStepValue(KM.NODAL_MAUX))
                ###node.SetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE, -1*node.GetSolutionStepValue(SMA.POINT_LOAD) / node.GetSolutionStepValue(KM.NODAL_MAUX))
                # print("For node id:",node.Id,", lagrange displacement=",node.GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT))
                # print("For node id:",node.Id,", for particle coupling force=",node.GetSolutionStepValue(VCA.DEMFEM_VOLUME_COUPLING_FORCE))
                #print("For node id:",node.Id,", point load=",node.GetSolutionStepValue(SMA.POINT_LOAD)) 
                # print("For node id:",node.Id,", velocity=",node.GetSolutionStepValue(KM.VELOCITY)) 
                
        # for node in self.model_part.Nodes: 
        #     print("For node id:",node.Id,", nodal coupling weight=",node.GetSolutionStepValue(VCA.NODAL_COUPLING_WEIGHT))
        #print(self.model_part.NumberOfConditions(0))

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : "",
            "penalty_max"           : 1e10,
            "timestep"              : 1e-5,
            "force_end_time"        : 1e-1,
            "boundary_tolerance"    : 1e-6,
            "y_fem_boundary"       : 0.04,
            "y_dem_boundary"       : 0.08,
            "weight_fem_boundary"  : 0.01,
            "weight_dem_boundary"  : 0.99
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults




