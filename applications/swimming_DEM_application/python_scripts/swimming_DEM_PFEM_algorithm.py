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
    
    def FluidInitialize(self):
        import kratos_pfem_fluid_ready_for_dem_coupling as solver_module
        self.solver_module = solver_module
        self.fluid_solver = self.solver_module.Solution()
        self.fluid_solver.Initialize()
        self.fluid_model_part = self.fluid_solver.main_model_part
        self.all_model_parts.Set('FluidPart', self.fluid_solver.main_model_part)
        
    def SetFluidSolverParameters(self):              
        self.pp = self.FluidSolverParameters()
        
    def CreateParts(self):
        # Order must be respected here
        # defining a fluid model
        self.all_model_parts.Add(ModelPart("FluidPart"))
        # defining a model part for the mixed part
        self.all_model_parts.Add(ModelPart("MixedPart"))
        
    class FluidSolverParameters(self):
        def __init__(self):
            pass
        


