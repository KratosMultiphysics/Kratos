from KratosMultiphysics import *
from KratosMultiphysics.TestPotentialApplication import *

import KratosMultiphysics

#import KratosMultiphysics.FluidDynamicsApplication
#import KratosMultiphysics.ConvectionDiffusionApplication 
#import KratosMultiphysics.FullPotentialFlowApplication
from array import array
from numpy import *


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyBCProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyBCProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "inlet_phi": 1.0,
                "velocity_infinity": [1.0,0.0,0],
                "density_infinity": 1.2,
                "skin_parts" : ["Inlet2D_Inlet_velocity_Auto1","Outlet2D_Outlet_pressure_Auto1","NoSlip2D_No_Slip_Auto1","NoSlip2D_No_Slip_Airfoil"]
            }  """ );
        
            
        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.velocity_infinity = [0,0,0]#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.density_infinity = settings["density_infinity"].GetDouble()
        self.inlet_phi = settings["inlet_phi"].GetDouble()
        
        #Find neighbours
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        neighbour_search = FindNodalNeighboursProcess(self.model_part,number_of_avg_elems,number_of_avg_nodes)
        neighbour_search.Execute()
        #create a new model part including all of the skin components
        self.all_skin_modelpart = self.model_part.CreateSubModelPart("all_skin_modelpart")
        #print(settings["skin_parts"].size())
        for part_id in range(settings["skin_parts"].size()):
            mp = Model[settings["skin_parts"][part_id].GetString()]
            for node in mp.Nodes:
              self.all_skin_modelpart.Nodes.append(node)
            for cond in mp.Conditions:
              self.all_skin_modelpart.Conditions.append(cond)
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.all_skin_modelpart, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        
        model_part = self.all_skin_modelpart
        
        #My stuff
        self.inlet_modelpart = self.model_part.CreateSubModelPart("inlet_modelpart")
        mp = Model[settings["skin_parts"][0].GetString()]
        #print('inlet_modelpart',mp)
        for node in mp.Nodes:
           self.inlet_modelpart.Nodes.append(node)
        for cond in mp.Conditions:           
           self.inlet_modelpart.Conditions.append(cond)
              
#              
        self.outlet_modelpart = self.model_part.CreateSubModelPart("outlet_modelpart")
        mp = Model[settings["skin_parts"][1].GetString()]
        for node in mp.Nodes:
           self.outlet_modelpart.Nodes.append(node)
        for cond in mp.Conditions:
           self.outlet_modelpart.Conditions.append(cond)
              
              
        self.slip_modelpart = self.model_part.CreateSubModelPart("slip_modelpart")
        mp = Model[settings["skin_parts"][2].GetString()]
        for node in mp.Nodes:
           self.slip_modelpart.Nodes.append(node)
        for cond in mp.Conditions:
           self.slip_modelpart.Conditions.append(cond)
           
           
        self.airfoil_modelpart = self.model_part.CreateSubModelPart("airfoil_modelpart")
        mp = Model[settings["skin_parts"][3].GetString()]
        for node in mp.Nodes:
           self.airfoil_modelpart.Nodes.append(node)
        for cond in mp.Conditions:
           self.airfoil_modelpart.Conditions.append(cond)

  
    def ExecuteInitializeSolutionStep(self):
         # DO I REALLY NEED THIS LINES?
         #for node in self.inlet_modelpart.Nodes:
         #   node.Free(KratosMultiphysics.TestFullPotentialApplication.POTENTIAL)   

         for cond in self.inlet_modelpart.Conditions:
           for node in cond.GetNodes():
             node.Fix(KratosMultiphysics.TestPotentialApplication.POTENTIAL)
             node.SetSolutionStepValue(KratosMultiphysics.TestPotentialApplication.POTENTIAL, 0, 2)
           
         for cond in self.outlet_modelpart.Conditions:
           cond.SetValue(TAU, self.velocity_infinity[0]*self.density_infinity)
         
         for cond in self.slip_modelpart.Conditions:
           cond.SetValue(TAU, 0)
           
         for cond in self.airfoil_modelpart.Conditions:
           cond.SetValue(TAU, 0)
