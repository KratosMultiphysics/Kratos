from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import swimming_DEM_algorithm 
import swimming_DEM_procedures as SDP

sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters
BaseAlgorithm = swimming_DEM_algorithm.Algorithm

class Algorithm(BaseAlgorithm):
                
    def SetFluidAlgorithm(self):
        import pfem_fluid_ready_for_dem_coupling as fluid_algorithm
        self.fluid_algorithm = fluid_algorithm.Solution()
        self.fluid_algorithm.main_path = self.main_path                
    
    def SetCouplingParameters(self, varying_parameters):         
     
        super(Algorithm,self).SetCouplingParameters(varying_parameters)       
        self.pp.domain_size = self.fluid_algorithm.ProjectParameters["problem_data"]["domain_size"].GetInt()

    def SetBetaParameters(self):

        self.pp.Dt = self.fluid_algorithm.GetDeltaTimeFromParameters()
        super(Algorithm,self).SetBetaParameters()

    #def SetCouplingParameters(self, varying_parameters):         
        #parameters_file = open("ProjectParametersDEM.json",'r')
        #self.pp.CFD_DEM = Parameters(parameters_file.read())
        #self.SetDoSolveDEMVariable()
        #self.pp.Dt = self.fluid_algorithm.GetDeltaTimeFromParameters()
        #self.SetBetaParameters()
        #self.SetCustomBetaParameters(varying_parameters)           
        #self.pp.domain_size = self.fluid_algorithm.ProjectParameters["problem_data"]["domain_size"].GetInt()
        #super(Algorithm,self).SetCouplingParameters(varying_parameters)
        
    def SetAllModelParts(self):
        self.all_model_parts = self.disperse_phase_algorithm.all_model_parts
        
        # defining a fluid model
        self.all_model_parts.Add(self.fluid_model_part, "FluidPart")        
        
        # defining a model part for the mixed part
        self.all_model_parts.Add(ModelPart("MixedPart"))  
        
        self.mixed_model_part = self.all_model_parts.Get('MixedPart')

    def Initialize(self):
        super(Algorithm,self).Initialize()
        self.TransferWallsFromPfemToDem()

    def TransferWallsFromPfemToDem(self):
        destination_model_part = self.disperse_phase_algorithm.rigid_face_model_part
        bodies_parts_list = self.fluid_algorithm.ProjectParameters["solver_settings"]["bodies_list"]
        for i in range(bodies_parts_list.size()):
            body_model_part_type = bodies_parts_list[i]["body_type"].GetString() 
            if body_model_part_type == "Rigid":
                body_model_part_name = bodies_parts_list[i]["body_name"].GetString()
                source_model_part = self.fluid_algorithm.main_model_part.GetSubModelPart(body_model_part_name) 
                SwimmingDemInPfemUtils().TransferWalls(source_model_part, destination_model_part)

    def FluidInitialize(self):
        
        self.fluid_algorithm.vars_man=self.vars_man        
        self.fluid_algorithm.Initialize() 
        bodies_parts_list = self.fluid_algorithm.ProjectParameters["solver_settings"]["bodies_list"]
        for i in range(bodies_parts_list.size()):
            body_model_part_type = bodies_parts_list[i]["body_type"].GetString() 
            if body_model_part_type == "Fluid":
                body_model_part_name = bodies_parts_list[i]["body_name"].GetString()
                self.fluid_model_part = self.fluid_algorithm.main_model_part.GetSubModelPart(body_model_part_name)
                break
        
    def TransferTimeToFluidSolver(self):
        if self.step < self.GetFirstStepForFluidComputation() or self.stationarity:
            self.fluid_algorithm.time = self.time
            #self.fluid_algorithm.step = self.step #DO NOT INCREASE STEP IN PFEM, IT CRASHES (PROBABLY IT MUST DO SPECIAL THINGS FOR STEP=0)
            #self.fluid_algorithm.main_model_part.ProcessInfo[STEP] = self.step #DO NOT INCREASE STEP IN PFEM, IT CRASHES (PROBABLY IT MUST DO SPECIAL THINGS FOR STEP=0)
            print(" [STEP:",self.step," TIME:",self.time,"]")
    
    def CloneTimeStep(self):
        
        if self.step < self.GetFirstStepForFluidComputation() or self.stationarity:
            self.fluid_algorithm.main_model_part.CloneTimeStep(self.time) 
        
    def FluidSolve(self, time = 'None'):
        
        self.fluid_algorithm.InitializeSolutionStep()
        self.fluid_algorithm.SolveSolutionStep()
        self.fluid_algorithm.FinalizeSolutionStep() 
        self.projection_module.UpdateDatabase(self.h_min)
        
    
    def GetFirstStepForFluidComputation(self):
        return 1;

    def SetCutsOutput(self):
        pass
    
    def SetDragOutput(self):
        pass
    
    def FinalizeDragOutput(self): 
        pass
    
    def SetPointGraphPrinter(self):
        pass
                
    def SetFluidSolverParameters(self):              
        
        self.pp.domain_size
        self.time           = self.pp.Start_time
        self.Dt             = self.pp.Dt
        self.out            = self.Dt
        Nsteps         = self.pp.nsteps
        if "REACTION" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.pp.fluid_vars)
        self.pp.variables_to_print_in_file
        if self.pp.type_of_inlet == 'ForceImposed':
            self.DEM_inlet = DEM_Force_Based_Inlet(self.DEM_inlet_model_part, self.pp.force)
          
    def SetPostUtils(self):
        general_model_part = self.fluid_algorithm.main_model_part
        self.post_utils = SDP.PostUtils(self.swimming_DEM_gid_io,
                                        self.pp,
                                        general_model_part,
                                        self.disperse_phase_algorithm.spheres_model_part,
                                        self.disperse_phase_algorithm.cluster_model_part,
                                        self.disperse_phase_algorithm.rigid_face_model_part,
                                        self.mixed_model_part)
    def SetEmbeddedTools(self):
        pass
            
    def PerformEmbeddedOperations(self):
        pass

        


