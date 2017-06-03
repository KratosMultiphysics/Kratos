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
        self.pp.Dt???????
        self.pp.problem_name,
        self.pp.VolumeOutput,
        self.pp.GiDPostMode,
        self.pp.GiDMultiFileFlag,
        self.pp.GiDWriteMeshFlag,
        self.pp.GiDWriteConditionsFlag
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
        
        
    class FluidSolverParameters(self):
        def __init__(self):
            pass
        


