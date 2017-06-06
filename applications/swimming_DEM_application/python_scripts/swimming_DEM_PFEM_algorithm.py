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
        
        self.pp.Dt = self.fluid_algorithm.GetDeltaTimeFromParameters()
        self.pp.domain_size = self.fluid_algorithm.ProjectParameters["problem_data"]["domain_size"].GetInt()
        super(Algorithm,self).SetCouplingParameters(varying_parameters)
        
    def SetAllModelParts(self):
        self.all_model_parts = self.disperse_phase_algorithm.all_model_parts
        
        # defining a fluid model
        self.all_model_parts.Add(ModelPart("FluidPart"))        
        
        #self.fluid_model_part = self.fluid_algorithm.main_model_part.GetSubModelPart("fluid_computing_domain")
        #self.fluid_model_part = self.fluid_algorithm.main_model_part.GetSubModelPart("Body1")
        self.all_model_parts.Set("FluidPart", self.fluid_model_part)
        
        # defining a model part for the mixed part
        self.all_model_parts.Add(ModelPart("MixedPart"))  
        
        self.mixed_model_part = self.all_model_parts.Get('MixedPart')
        
    def FluidInitialize(self):
        
        self.fluid_algorithm.vars_man=self.vars_man
        self.vars_man.AddExtraProcessInfoVariablesToFluidModelPart(self.pp, self.fluid_model_part)
        self.fluid_algorithm.Initialize()   
        self.fluid_model_part = self.fluid_algorithm.main_model_part.GetSubModelPart("Body1")
        
    def CloneTimeStep(self):
        pass
        
    def FluidSolve(self, time = 'None'):
        
        self.fluid_algorithm.InitializeSolutionStep()
        self.fluid_algorithm.SolveSolutionStep()
        self.fluid_algorithm.FinalizeSolutionStep()                
        
    def SetCutsOutput(self):
        pass
    
    def SetDragOutput(self):
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

        


