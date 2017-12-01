from KratosMultiphysics import *
import swimming_DEM_procedures as SDP
import swimming_DEM_algorithm
BaseAlgorithm = swimming_DEM_algorithm.Algorithm
import math
import chandelier as ch
import chandelier_parameters as ch_pp

def SafeCrossProduct(a, b):
    c0 = a[1]*b[2] - a[2]*b[1]
    c1 = a[2]*b[0] - a[0]*b[2]
    c2 = a[0]*b[1] - a[1]*b[0]
    return c0, c1, c2

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
        is_rotating_frame = self.pp.CFD_DEM["frame_of_reference_type"].GetInt()
        omega = self.pp.CFD_DEM["angular_velocity_of_frame_Z"].GetDouble()

        for node in self.spheres_model_part.Nodes:
            coor_calculated = [node.X, node.Y, node.Z]

            if is_rotating_frame:
                coor_calculated[0] = node.X * math.cos(omega * time) - node.Y * math.sin(omega * time)
                coor_calculated[1] = node.X * math.sin(omega * time) + node.Y * math.cos(omega * time)

            self.radial_error = self.results_database.CalculateError(time, coor_calculated)
            self.error_time = time

    def PerformZeroStepInitializations(self):
        # Impose initial velocity to be the terminal velocity
        ch.sim.CalculateNonDimensionalVars()
        terminal_velocity = ch.sim.NDw0 *  ch_pp.R *  ch_pp.omega

        for node in self.spheres_model_part.Nodes:
            node_to_follow_id = node.Id
            node.SetSolutionStepValue(VELOCITY_Z, terminal_velocity)

        terminal_velocity_z = 2. / 9 * 9.81 * ch_pp.a ** 2 / (ch_pp.nu * ch_pp.rho_f) * (ch_pp.rho_f - ch_pp.rho_p)

        for node in self.spheres_model_part.Nodes:
            r = [node.X, node.Y, node.Z]
            cross_omega_r = [0., 0., 0.]

            if self.pp.CFD_DEM["frame_of_reference_type"].GetInt():
                omega_frame = [0, 0, self.pp.CFD_DEM["angular_velocity_of_frame_Z"].GetDouble()]
                cross_omega_r = list(SafeCrossProduct(omega_frame, r))

            node.SetSolutionStepValue(VELOCITY_X, ch_pp.u0 - cross_omega_r[0])
            node.SetSolutionStepValue(VELOCITY_Y, ch_pp.v0 - cross_omega_r[1])
            node.SetSolutionStepValue(VELOCITY_Z, terminal_velocity_z)
            node.Fix(VELOCITY_Z)
            node.SetSolutionStepValue(VELOCITY_OLD_X, ch_pp.u0 - cross_omega_r[0])
            node.SetSolutionStepValue(VELOCITY_OLD_Y, ch_pp.v0 - cross_omega_r[1])
            node.SetSolutionStepValue(VELOCITY_OLD_Z, terminal_velocity_z)
            node.Fix(VELOCITY_OLD_Z)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, ch_pp.u0 - cross_omega_r[0])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, ch_pp.v0 - cross_omega_r[1])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)

    def ApplyForwardCoupling(self, alpha = 'None'):
        self.projection_module.ApplyForwardCoupling(alpha)
        is_rotating_frame = self.pp.CFD_DEM["frame_of_reference_type"].GetInt()
        for node in self.spheres_model_part.Nodes:
            omega = ch_pp.omega
            r = [node.X, node.Y, node.Z]
            vx = - omega * r[1]
            vy =   omega * r[0]
            ax = - r[0] * omega ** 2
            ay = - r[1] * omega ** 2
            v = [vx, vy, 0]
            cross_omega_r = [0., 0., 0.]
            cross_omega_v = [0., 0., 0.]
            cross_omega_omega_r = [0., 0., 0.]

            if is_rotating_frame:
                omega_frame = [0, 0, self.pp.CFD_DEM["angular_velocity_of_frame_Z"].GetDouble()]
                cross_omega_r = list(SafeCrossProduct(omega_frame, r))
                cross_omega_omega_r = list(SafeCrossProduct(omega_frame, cross_omega_r))
                v = [v[i] - cross_omega_r[i] for i in range(3)]
                cross_omega_v = list(SafeCrossProduct(omega_frame, v))

            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, v[0])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, v[1])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_X, ax - 2 * cross_omega_v[0] - cross_omega_omega_r[0])
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Y, ay - 2 * cross_omega_v[1] - cross_omega_omega_r[1])
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Z, 0.0)

    def ApplyForwardCouplingOfVelocityOnly(self, time):
        is_rotating_frame = self.pp.CFD_DEM["frame_of_reference_type"].GetInt()

        for node in self.spheres_model_part.Nodes:
            self.projection_module.ApplyForwardCouplingOfVelocityOnly()
            r = [node.X, node.Y, node.Z]
            self.results_database.MakeReading(time, r)
            new_vx = - ch_pp.omega * r[1]
            new_vy =   ch_pp.omega * r[0]

            cross_omega_r = [0., 0., 0.]

            if is_rotating_frame:
                omega_frame = [0, 0, self.pp.CFD_DEM["angular_velocity_of_frame_Z"].GetDouble()]
                cross_omega_r = SafeCrossProduct(omega_frame, r)

            node.SetSolutionStepValue(SLIP_VELOCITY_X, new_vx - cross_omega_r[0])
            node.SetSolutionStepValue(SLIP_VELOCITY_Y, new_vy - cross_omega_r[1])

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
