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
    def __init__(self):
        import time as timer
        self.timer = timer
        self.simulation_start_time = timer.time()
        BaseAlgorithm.__init__(self)

    def Initialize(self):
        BaseAlgorithm.Initialize(self)
        self.fluid_solver.Initialize()
        print("fluid solver created")
        sys.stdout.flush()

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
