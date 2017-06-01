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
    
    def ReadModelParts(self, starting_node_Id = 0, starting_elem_Id = 0, starting_cond_Id = 0):
        self.ReadFluidModelPart()
        fluid_mp = self.all_model_parts.Get('FluidPart')
        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(fluid_mp)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(fluid_mp)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(fluid_mp)
        self.ReadDEMModelParts(max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)

    def ReadFluidModelPart(self):
        pass        

    def Initialize(self):
        BaseAlgorithm.Initialize(self)

    def FluidInitialize(self):
        self.fluid_solver.Initialize()

    def AddExtraVariables(self):
        spheres_model_part = self.all_model_parts.Get('SpheresPart')
        fluid_model_part = self.all_model_parts.Get('FluidPart')

                # building lists of variables for which memory is to be allocated
        # TEMPORARY, HORRIBLE !!!
        import variables_management as vars_man

        vars_man.ConstructListsOfVariables(self.pp)
        #_____________________________________________________________________________________________________________________________________
        #
        #                               F L U I D    B L O C K    B E G I N S
        #_____________________________________________________________________________________________________________________________________

        # defining variables to be used
        # GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
        # CODE IS NOT USING IT

        fluid_model_part = self.all_model_parts.Get('FluidPart')

        if "REACTION" in self.pp.nodal_results:
            fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in self.pp.nodal_results:
            fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

        # importing the solvers needed
        SolverSettings = self.pp.FluidSolverConfiguration
        self.solver_module = import_solver(SolverSettings)

        # caution with breaking up this block (memory allocation)! {
        self.solver_module.AddVariables(fluid_model_part, SolverSettings)
        vars_man.AddNodalVariables(fluid_model_part, self.pp.fluid_vars)  #     MOD.
        # }
        # self.ReadFluidModelPart()
        # Creating necessary directories
        [post_path, data_and_results, graphs_path, MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(self.pp.CFD_DEM.problem_name))

        #_____________________________________________________________________________________________________________________________________
        #
        #                               F L U I D    B L O C K    E N D S
        #_____________________________________________________________________________________________________________________________________

        # Add variables

        vars_man.AddNodalVariables(spheres_model_part, self.pp.dem_vars)
        vars_man.AddNodalVariables(self.rigid_face_model_part, self.pp.rigid_faces_vars)
        vars_man.AddNodalVariables(self.DEM_inlet_model_part, self.pp.inlet_vars)
        # adding extra process info variables
        vars_man.AddingExtraProcessInfoVariables(self.pp, fluid_model_part, spheres_model_part)

    def SetFluidBufferSizeAndAddAdditionalDofs(self):
        spheres_model_part = self.all_model_parts.Get('SpheresPart')
        fluid_model_part = self.all_model_parts.Get('FluidPart')
        SolverSettings = self.pp.FluidSolverConfiguration

        fluid_model_part.SetBufferSize(3)
        self.solver_module.AddDofs(fluid_model_part, SolverSettings)
        SDP.AddExtraDofs(self.pp, fluid_model_part, spheres_model_part, self.cluster_model_part, self.DEM_inlet_model_part)
        os.chdir(self.main_path)

        self.KRATOSprint("\nInitializing Problem...")

    

    def FluidSolve(self, time = 'None'):
        self.fluid_solver.Solve()


    def FillHistoryForcePrecalculatedVectors(self):
        # Warning: this estimation is based on a constant time step for DEM. This is usually the case, but could not be so. A more robust implementation is needed!
        N_steps = int(self.pp.CFD_DEM.FinalTime / self.pp.CFD_DEM.MaxTimeStep) + 20
        spheres_model_part = self.all_model_parts.Get('SpheresPart')
        if self.pp.CFD_DEM.basset_force_type > 0:
            self.basset_force_tool.FillDaitcheVectors(N_steps, self.pp.CFD_DEM.quadrature_order, self.pp.CFD_DEM.time_steps_per_quadrature_step)
        if self.pp.CFD_DEM.basset_force_type >= 3 or self.pp.CFD_DEM.basset_force_type == 1:
            self.basset_force_tool.FillHinsbergVectors(spheres_model_part, self.pp.CFD_DEM.number_of_exponentials, self.pp.CFD_DEM.number_of_quadrature_steps_in_window)

