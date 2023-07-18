# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA 

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

import math

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
        self.tolerance= self.settings["velocity_tolerance"].GetDouble() 
        self.y_fem_boundary= self.settings["y_fem_boundary"].GetDouble() 
        self.y_dem_boundary= self.settings["y_dem_boundary"].GetDouble() 
        self.weight_fem_boundary = 0.9
        self.weight_dem_boundary = 0.1
        self.step=0
        self.timesteps=self.force_end_time/self.dt
        self.max_force=0.025
        self.force_slope=self.max_force/self.timesteps # to reach max force within force_end_time


    def InitializeCouplingIteration(self):


        self.step+=1
        pointload=min(self.force_slope*self.step,self.max_force) # gradually increasing point load

        node_ids = [1,2,3,4]  # for assigning point loads to the top nodes

        for node in self.model_part.Nodes:  
          
                if node.Id in node_ids:   # assigning point loads

                    node.SetSolutionStepValue(SMA.POINT_LOAD,[0,-pointload,0])
                    print("For node id:",node.Id,", point load=",node.GetSolutionStepValue(SMA.POINT_LOAD))

                if node.Y > self.y_fem_boundary and node.Y < self.y_dem_boundary: # assigning weights to the nodes in the hybrid region
                    weight= self.weight_fem_boundary + (self.weight_dem_boundary-self.weight_fem_boundary)*(node.Y-self.y_fem_boundary)/(self.y_dem_boundary-self.y_fem_boundary)
                    node.SetSolutionStepValue(SMA.NODAL_COUPLING_WEIGHT, weight)
               
                if  node.Y<self.y_fem_boundary and node.Y>self.y_dem_boundary: # assigning weights to the nodes in the hybrid region
                    weight= self.weight_dem_boundary + (self.weight_fem_boundary-self.weight_dem_boundary)*(node.Y-self.y_dem_boundary)/(self.y_fem_boundary-self.y_dem_boundary)
                    node.SetSolutionStepValue(SMA.NODAL_COUPLING_WEIGHT, weight)


                total_mass = node.GetSolutionStepValue(KM.NODAL_MAUX)

                if(total_mass!=0): # check for hybrid region
                    node.SetSolutionStepValue(SMA.POINT_LOAD,[0,0,0]) #nodal coupling forces set to zero before every iteration
                    node.SetSolutionStepValue(KM.VELOCITY_LAPLACIAN,[0,0,0]) # force to map to dem (set to zero before every iteration)
                    Displacement_dem = node.GetSolutionStepValue(KM.LINEAR_MOMENTUM) / total_mass # total lagrange method- calculating homogenised displacement
                    node.SetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT, Displacement_dem - node.GetSolutionStepValue(KM.DISPLACEMENT))# total lagrange method : calcualting displacement difference

        for elem in self.model_part.Elements:
             if(elem.GetNodes()[0].GetSolutionStepValue(KM.NODAL_MAUX))!=0:  
                for i in range(elem.GetGeometry().IntegrationPointsNumber()): #gauss quadrature
                    w = 1 
                    J = (elem.GetGeometry().DeterminantOfJacobian(i))
                    shape_functions = elem.GetGeometry().ShapeFunctionsValues()
                    for n in range(len(elem.GetNodes())):
                        for m in range(len(elem.GetNodes())):
                            vol = self.penalty_max * w * J * shape_functions[(i,n)] * shape_functions[(i,m)]
                            elem.GetNodes()[n].SetSolutionStepValue(SMA.POINT_LOAD, elem.GetNodes()[n].GetSolutionStepValue(SMA.POINT_LOAD) + 2* vol * elem.GetNodes()[n].GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT))
                            if(elem.GetNodes()[n].GetSolutionStepValue(KM.NODAL_MAUX))>0.002: # mass of one particle is 0.00131kg, checking if the node has contribution from 2 elements
                                 elem.GetNodes()[n].SetSolutionStepValue(KM.VELOCITY_LAPLACIAN,elem.GetNodes()[n].GetSolutionStepValue(KM.VELOCITY_LAPLACIAN) - vol * elem.GetNodes()[n].GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT)) #storing force to map to dem 
                            else:
                                 elem.GetNodes()[n].SetSolutionStepValue(KM.VELOCITY_LAPLACIAN,elem.GetNodes()[n].GetSolutionStepValue(KM.VELOCITY_LAPLACIAN) - 2 * vol * elem.GetNodes()[n].GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT)) #storing force to map to dem   
     
        print("After calculation of point loads")
        for node in self.model_part.Nodes: 
                print("For node id:",node.Id,", point load=",node.GetSolutionStepValue(SMA.POINT_LOAD)) 

        print(self.model_part.NumberOfConditions(0))

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : "",
            "penalty_max"           : 1e11,
            "timestep"              : 1e-5,
            "force_end_time"        : 1e-1,
            "velocity_tolerance"    : 1e-6
            "y_fem_boundary"       : 0.16
            "y_dem_boundary"       : 0.08
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults




