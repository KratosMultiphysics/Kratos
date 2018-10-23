from KratosMultiphysics import *
import swimming_DEM_procedures as SDP
import swimming_DEM_algorithm
BaseAlgorithm = swimming_DEM_algorithm.Algorithm
import math
import chandelier as ch
import chandelier_parameters as ch_pp

def Cross(a, b):
    c0 = a[1]*b[2] - a[2]*b[1]
    c1 = a[2]*b[0] - a[0]*b[2]
    c2 = a[0]*b[1] - a[1]*b[0]
    return Vector([c0, c1, c2])


class Algorithm(BaseAlgorithm):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, model, varying_parameters)

    def GetVelocityRelativeToMovingFrame(self, r_rel, v_glob):
        cross_omega_r = Cross(self.frame_angular_vel, r_rel)
        return v_glob - cross_omega_r

    def GetAccelerationRelativeToMovingFrame(self, r_rel, v_rel, a_glob):
        cross_omega_r = Cross(self.frame_angular_vel, r_rel)
        cross_omega_v = Cross(self.frame_angular_vel, v_rel)
        cross_omega_omega_r = Cross(self.frame_angular_vel, cross_omega_r)
        return a_glob - 2 * cross_omega_v - cross_omega_omega_r

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
        self.frame_angular_vel = Vector([0, 0, self.pp.CFD_DEM["angular_velocity_of_frame_Z"].GetDouble()])
        self.is_rotating_frame = self.pp.CFD_DEM["frame_of_reference_type"].GetInt()

    def SetUpResultsDatabase(self):
        import candelier_hdf5
        self.results_database = candelier_hdf5.ResultsCandelier(self.pp, self.main_path)

    def DEMSolve(self, time = 'None'):
        self.disperse_phase_solution.solver.Solve()
        omega = self.pp.CFD_DEM["angular_velocity_of_frame_Z"].GetDouble()

        for node in self.spheres_model_part.Nodes:
            coor_calculated = [node.X, node.Y, node.Z]

            if self.is_rotating_frame:
                sin = math.sin(omega * time)
                cos = math.cos(omega * time)
                coor_calculated[0] = node.X * cos - node.Y * sin
                coor_calculated[1] = node.X * sin + node.Y * cos

            self.radial_error = self.results_database.CalculateError(time, coor_calculated)
            self.error_time = time

    def PerformZeroStepInitializations(self):
        # Impose initial velocity to be that of the fluid for the x/y-components
        # and the terminal velocity for the z-component
        ch.sim.CalculateNonDimensionalVars()
        terminal_velocity_z = 2. / 9 * 9.81 * ch_pp.a ** 2 / (ch_pp.nu * ch_pp.rho_f) * (ch_pp.rho_f - ch_pp.rho_p)

        for node in self.spheres_model_part.Nodes:
            r = Vector([node.X, node.Y, node.Z])
            v0 = Vector([ch_pp.u0, ch_pp.v0, terminal_velocity_z])

            if self.is_rotating_frame:
                v0 = self.GetVelocityRelativeToMovingFrame(r_rel = r, v_glob = v0)

            node.SetSolutionStepValue(VELOCITY, v0)
            node.Fix(VELOCITY_Z)
            node.SetSolutionStepValue(VELOCITY_OLD, v0)
            node.Fix(VELOCITY_OLD_Z)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, v0[0])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, v0[1])
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)

    def ApplyForwardCoupling(self, alpha = 'None'):
        self.projection_module.ApplyForwardCoupling(alpha)

        for node in self.spheres_model_part.Nodes:
            omega = ch_pp.omega
            r = Vector([node.X, node.Y, node.Z])
            vx = - omega * r[1]
            vy =   omega * r[0]
            v = Vector([vx, vy, 0])
            ax = - r[0] * omega ** 2
            ay = - r[1] * omega ** 2
            a = Vector([ax, ay, 0])

            if self.is_rotating_frame:
                v = self.GetVelocityRelativeToMovingFrame(r_rel = r, v_glob = v)
                a = self.GetAccelerationRelativeToMovingFrame(r_rel = r, v_rel = v, a_glob = a)

            node.SetSolutionStepValue(FLUID_VEL_PROJECTED, v)
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED, a)

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self, time):
        for node in self.spheres_model_part.Nodes:
            r = Vector([node.X, node.Y, node.Z])
            self.results_database.MakeReading(time, r)
            new_vx = - ch_pp.omega * r[1]
            new_vy =   ch_pp.omega * r[0]
            new_v = Vector([new_vx, new_vy, 0.])

            if self.is_rotating_frame:
                new_v = self.GetVelocityRelativeToMovingFrame(r_rel = r, v_glob = new_v)

            # the current FLUID_VEL_PROJECTED is still needed and so we use
            # SLIP_VELOCITY to store it.
            node.SetSolutionStepValue(SLIP_VELOCITY, new_v)

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
