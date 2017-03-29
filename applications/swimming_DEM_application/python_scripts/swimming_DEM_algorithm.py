from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import main_script as DEM_algorithm
import swimming_DEM_procedures as SDP

sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters
BaseAlgorithm = DEM_algorithm.Solution

class Algorithm(BaseAlgorithm):
    def __init__(self, pp):
        self.pp = pp
        import time as timer
        self.timer = timer
        self.simulation_start_time = timer.time()
        BaseAlgorithm.__init__(self)
        # defining a model part for the mixed part
        self.mixed_model_part = ModelPart("MixedPart")

    def Initialize(self):
        BaseAlgorithm.Initialize(self)
        self.fluid_solver.Initialize()
        os.chdir(self.pp.main_path)
        model_part_io_fluid = ModelPartIO(self.pp.problem_name)
        model_part_io_fluid.ReadModelPart(self.all_model_parts.Get('FluidPart'))
        print("fluid solver created")
        sys.stdout.flush()

    def AddExtraVariables(self):

        import variables_management as vars_man
        fluid_model_part = self.all_model_parts.Get('FluidPart')

                # building lists of variables for which memory is to be allocated
        # TEMPORARY, HORRIBLE !!!
        vars_man.ConstructListsOfVariables(self.pp)
        #_____________________________________________________________________________________________________________________________________
        #
        #                               F L U I D    B L O C K    B E G I N S
        #_____________________________________________________________________________________________________________________________________

        # defining variables to be used
        # GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
        # CODE IS NOT USING IT

        variables_dictionary = {"PRESSURE"   : PRESSURE,
                                "VELOCITY"   : VELOCITY,
                                "MU"         : MU,         #    MOD.
                                "BUOYANCY"   : BUOYANCY,   #    MOD.
                                "DRAG_FORCE" : DRAG_FORCE,  #    MOD.
                                "LIFT_FORCE" : LIFT_FORCE} #    MOD.

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



        # Creating necessary directories
        [post_path, data_and_results, graphs_path, MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(self.pp.CFD_DEM.problem_name))

        #_____________________________________________________________________________________________________________________________________
        #
        #                               F L U I D    B L O C K    E N D S
        #_____________________________________________________________________________________________________________________________________

        # Add variables

        vars_man.AddNodalVariables(self.spheres_model_part, self.pp.dem_vars)
        vars_man.AddNodalVariables(self.rigid_face_model_part, self.pp.rigid_faces_vars)
        vars_man.AddNodalVariables(self.DEM_inlet_model_part, self.pp.inlet_vars)
        # adding extra process info variables
        vars_man.AddingExtraProcessInfoVariables(self.pp, fluid_model_part, self.spheres_model_part)

        fluid_model_part.SetBufferSize(3)
        self.solver_module.AddDofs(fluid_model_part, SolverSettings)
        SDP.AddExtraDofs(self.pp, fluid_model_part, self.spheres_model_part, self.cluster_model_part, self.DEM_inlet_model_part)

        os.chdir(self.main_path)

        self.KRATOSprint("\nInitializing Problem...")


    def SetSolverStrategy(self):
        import swimming_sphere_strategy as SolverStrategy
        return SolverStrategy

    def SelectScheme(self):
        scheme = BaseAlgorithm.SelectScheme(self)
        if DEM_parameters.IntegrationScheme == 'Hybrid_Bashforth':
            return HybridBashforthScheme()
        elif DEM_parameters.ElementType == "SwimmingNanoParticle":
            return TerminalVelocityScheme()
        else:
            return None

    def SetSolver(self):
        return self.solver_strategy.SwimmingStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, DEM_parameters, self.procedures)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        fluid_mp = self.all_model_parts.Get('FluidPart')
        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(fluid_mp)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(fluid_mp)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(fluid_mp)
        BaseAlgorithm.ReadModelParts(self, max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)

    def TellTime(self, time):
        print("\n", "TIME = ", time)
        print('ELAPSED TIME = ', self.timer.time() - self.simulation_start_time)
        print()
        sys.stdout.flush()

    def TellFinalSummary(self, time, step, DEM_step):
        print()
        print()
        print('*************************************************************')
        print("CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.")
        simulation_elapsed_time = self.timer.time() - self.simulation_start_time
        print("Elapsed time: " + "%.5f"%(simulation_elapsed_time) + " s ")
        print("per fluid time step: " + "%.5f"%(simulation_elapsed_time / step) + " s ")
        print("per DEM time step: " + "%.5f"%(simulation_elapsed_time / DEM_step) + " s")
        print('*************************************************************')
        print()
        sys.stdout.flush()

    def GetFluidSolveCounter(self, pp):
        return SDP.Counter(is_dead = (pp.CFD_DEM.drag_force_type == 9))

    def GetEmbeddedCounter(self, pp):
        return SDP.Counter(1, 3, pp.CFD_DEM.embedded_option)  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self, pp):
        return SDP.Counter(1, 1, pp.CFD_DEM.coupling_level_type > 1)

    def GetBackwardCouplingCounter(self, pp):
        return SDP.Counter(1, 1, pp.CFD_DEM.coupling_level_type or self.pp.CFD_DEM.print_PRESSURE_GRADIENT_option)

    def GetStationarityCounter(self, pp):
        return SDP.Counter(pp.CFD_DEM.time_steps_per_stationarity_step, 1, pp.CFD_DEM.stationary_problem_option)

    def GetPrintCounter(self, pp):
        return SDP.Counter(1, 1, 10) # still unused

    def GetDebugInfo(self, pp):
        return SDP.Counter(pp.CFD_DEM.debug_tool_cycle, 1, pp.CFD_DEM.print_debug_info_option)

    def GetParticlesResultsCounter(self, pp):
        return SDP.Counter(pp.CFD_DEM.print_particles_results_cycle, 1, pp.CFD_DEM.print_particles_results_option)

    def HistoryForceQuadratureCounter(self, pp):
        return SDP.Counter(pp.CFD_DEM.time_steps_per_quadrature_step, 1, pp.CFD_DEM.basset_force_type)

    def GetRunCode(self, pp):
        return SDP.CreateRunCode(pp)

    def ActivateTurbulenceModel(self):
        fluid_model_part = self.all_model_parts.Get('FluidPart')

        if self.pp.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras":
            # apply the initial turbulent viscosity on all of the nodes
            turb_visc = self.pp.FluidSolverConfiguration.TurbulentViscosity
            for node in fluid_model_part.Nodes:
                node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
                visc = node.GetSolutionStepValue(VISCOSITY)
                node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
                if node.IsFixed(VELOCITY_X) and node.GetSolutionStepValue(VELOCITY_X, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Y) and node.GetSolutionStepValue(VELOCITY_Y, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Z) and node.GetSolutionStepValue(VELOCITY_Z, 0) != 0.0:
                    node.Fix(TURBULENT_VISCOSITY)

            # select nodes on the wall
            self.fluid_solver.wall_nodes = []
            for i in SolverSettings.SA_wall_group_ids:
                # get the nodes of the wall for SA.
                nodes = fluid_model_part.GetNodes(i)
                for node in nodes:
                    self.fluid_solver.wall_nodes.append(node)
                    node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
                    node.Fix(TURBULENT_VISCOSITY)
