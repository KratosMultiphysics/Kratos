from KratosMultiphysics import *
import swimming_DEM_procedures as SDP
import swimming_DEM_algorithm
BaseAlgorithm = swimming_DEM_algorithm.Algorithm
import math
import chandelier as ch
import chandelier_parameters as ch_pp

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead = True)

    def GetEmbeddedCounter(self):
        return SDP.Counter(1, 3, self.pp.CFD_DEM["embedded_option"].GetBool())  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 4, 0)

    def GetDebugInfo(self):
        return SDP.Counter(self.pp.CFD_DEM["debug_tool_cycle"].GetInt(), 1, is_dead = 1)

    def SetCustomBetaParameters(self, custom_parameters): # These are input parameters that have not yet been transferred to the interface
        BaseAlgorithm.SetCustomBetaParameters(self, custom_parameters)
        ch_pp.include_history_force = bool(self.pp.CFD_DEM["basset_force_type"].GetInt())
        ch.sim = ch.AnalyticSimulator(ch_pp)

    def SetUpResultsDatabase(self):
        import candelier_hdf5
        self.results_database = candelier_hdf5.ResultsCandelier(self.pp, self.main_path)

    def DEMSolve(self, time = 'None'):
        self.disperse_phase_solution.solver.Solve()
        for node in self.spheres_model_part.Nodes:
            coor_calculated = [node.X, node.Y, node.Z]
            self.radial_error = self.results_database.CalculateError(time, coor_calculated)
            self.error_time = time

    def PerformZeroStepInitializations(self):
        # Impose initial velocity to be the terminal velocity
        ch.sim.CalculateNonDimensionalVars()
        terminal_velocity = ch.sim.NDw0 *  ch_pp.R *  ch_pp.omega

        for node in self.spheres_model_part.Nodes:
            node_to_follow_id = node.Id
            node.SetSolutionStepValue(VELOCITY_Z, terminal_velocity)

        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, ch_pp.u0)
            node.SetSolutionStepValue(VELOCITY_Y, ch_pp.v0)
            node.SetSolutionStepValue(VELOCITY_Z, 2. / 9 * 9.81 * ch_pp.a ** 2 / (ch_pp.nu * ch_pp.rho_f) * (ch_pp.rho_f - ch_pp.rho_p))
            node.Fix(VELOCITY_Z)
            node.SetSolutionStepValue(VELOCITY_OLD_X, ch_pp.u0)
            node.SetSolutionStepValue(VELOCITY_OLD_Y, ch_pp.v0)
            node.SetSolutionStepValue(VELOCITY_OLD_Z, 2. / 9 * 9.81 * ch_pp.a ** 2 / (ch_pp.nu * ch_pp.rho_f) * (ch_pp.rho_f - ch_pp.rho_p))
            node.Fix(VELOCITY_OLD_Z)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, ch_pp.u0)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, ch_pp.v0)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)

    def ApplyForwardCoupling(self, alpha = 'None'):
        self.projection_module.ApplyForwardCoupling(alpha)
        for node in self.spheres_model_part.Nodes:
            x = node.X
            y = node.Y
            z = node.Z
            r = math.sqrt(x ** 2 + y ** 2)
            omega = ch_pp.omega
            vx = - omega * y
            vy =   omega * x
            ax = - x * omega ** 2
            ay = - y * omega ** 2
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, vx)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, vy)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_X, ax)
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Y, ay)
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Z, 0.0)

    def ApplyForwardCouplingOfVelocityOnly(self, time):
        for node in self.spheres_model_part.Nodes:
            self.projection_module.ApplyForwardCouplingOfVelocityOnly()
            self.x = node.X
            self.y = node.Y
            self.z = node.Z
            self.results_database.MakeReading(time, [self.x, self.y, self.z])
            r = math.sqrt(self.x ** 2 + self.y ** 2)
            new_vx = - ch_pp.omega * self.y
            new_vy =   ch_pp.omega * self.x
            node.SetSolutionStepValue(SLIP_VELOCITY_X, new_vx)
            node.SetSolutionStepValue(SLIP_VELOCITY_Y, new_vy)

    def PerformFinalOperations(self, time = None):
        self.results_database.WriteToHDF5()
        dt_quad_over_dt = self.pp.CFD_DEM["delta_time_quadrature"].GetDouble() / self.pp.CFD_DEM["MaxTimeStep"].GetDouble()
        os.chdir(self.main_path)
        sys.stdout.path_to_console_out_file
        import shutil
        folder_name = self.post_path + '_FINISHED_AT_t=' + str(round(time, 1))
        try:
            shutil.rmtree(folder_name)
        except OSError:
            pass

    def GetReturnValue(self):
        return self.radial_error
