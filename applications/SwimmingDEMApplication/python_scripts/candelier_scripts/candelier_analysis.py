from KratosMultiphysics import *
import swimming_DEM_procedures as SDP
import math
import candelier_scripts.candelier as candelier
import candelier_scripts.candelier_parameters as candelier_pp

def Cross(a, b):
    c0 = a[1]*b[2] - a[2]*b[1]
    c1 = a[2]*b[0] - a[0]*b[2]
    c2 = a[0]*b[1] - a[1]*b[0]
    return Vector([c0, c1, c2])

from swimming_DEM_analysis import SwimmingDEMAnalysis

class CandelierBenchmarkAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super(CandelierBenchmarkAnalysis, self).__init__(model, varying_parameters)
        self._GetSolver().is_rotating_frame = self.pp.CFD_DEM["frame_of_reference_type"].GetInt()
        self.disperse_phase_solution.mdpas_folder_path = os.path.join(self.disperse_phase_solution.main_path, 'candelier_tests')

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead = True)

    def GetEmbeddedCounter(self):
        return SDP.Counter(1, 3, self.pp.CFD_DEM["embedded_option"].GetBool())  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 4, 0)

    def GetDebugInfo(self):
        return SDP.Counter(self.pp.CFD_DEM["debug_tool_cycle"].GetInt(), 1, is_dead = 1)

    def SetCustomBetaParameters(self, custom_parameters): # These are input parameters that have not yet been transferred to the interface
        super(CandelierBenchmarkAnalysis, self).SetCustomBetaParameters(custom_parameters)
        candelier_pp.include_history_force = bool(self.pp.CFD_DEM["basset_force_type"].GetInt())
        candelier.sim = candelier.AnalyticSimulator(candelier_pp)
        self.pp.CFD_DEM["fluid_already_calculated"].SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("load_derivatives").SetBool(False)

    def SetUpResultsDatabase(self):
        import candelier_scripts.candelier_hdf5 as candelier_hdf5
        self.results_database = candelier_hdf5.ResultsCandelier(self.pp, self.main_path)

    def PerformZeroStepInitializations(self):
        # Impose initial velocity to be that of the fluid for the x/y-components
        # and the terminal velocity for the z-component
        candelier.sim.CalculateNonDimensionalVars()
        terminal_velocity_z = 2. / 9 * 9.81 * candelier_pp.a ** 2 / (candelier_pp.nu * candelier_pp.rho_f) * (candelier_pp.rho_f - candelier_pp.rho_p)

        for node in self.spheres_model_part.Nodes:
            r = Vector([node.X, node.Y, node.Z])
            v0 = Vector([candelier_pp.u0, candelier_pp.v0, terminal_velocity_z])

            if self._GetSolver().is_rotating_frame:
                v0 = self._GetSolver().GetVelocityRelativeToMovingFrame(r_rel = r, v_glob = v0)

            node.SetSolutionStepValue(VELOCITY, v0)
            node.Fix(VELOCITY_Z)
            node.SetSolutionStepValue(VELOCITY_OLD, v0)
            node.Fix(VELOCITY_OLD_Z)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, v0[0])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, v0[1])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, v0[2])

    def _CreateSolver(self):
        import candelier_scripts.candelier_dem_solver as sdem_solver
        self.pp.field_utility = self.GetFieldUtility()
        return sdem_solver.CandelierDEMSolver(self.model,
                                              self.pp.CFD_DEM,
                                              self.fluid_solution._GetSolver(),
                                              self.disperse_phase_solution._GetSolver(),
                                              self.pp)

    def FinalizeSolutionStep(self):
        super(CandelierBenchmarkAnalysis, self).FinalizeSolutionStep()
        if self.pp.CFD_DEM["do_print_results_option"].GetBool():
            for node in self.spheres_model_part.Nodes:
                r = Vector([node.X, node.Y, node.Z])
                self.results_database.MakeReading(self.time, r)

                coor_calculated = [node.X, node.Y, node.Z]

                if self._GetSolver().is_rotating_frame:
                    sin = math.sin(self._GetSolver().omega * self.time)
                    cos = math.cos(self._GetSolver().omega * self.time)
                    coor_calculated[0] = node.X * cos - node.Y * sin
                    coor_calculated[1] = node.X * sin + node.Y * cos

                self.radial_error = self.results_database.CalculateError(self.time, coor_calculated)
                self.error_time = self.time