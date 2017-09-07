from KratosMultiphysics import *
import swimming_DEM_procedures as SDP
import swimming_DEM_algorithm
BaseAlgorithm = swimming_DEM_algorithm.Algorithm
import math

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = dict()):
        BaseAlgorithm.__init__(self, varying_parameters)

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead = True)

    def GetEmbeddedCounter(self):
        return SDP.Counter(1, 3, self.pp.CFD_DEM["embedded_option"].GetBool())  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 4, 0)

    def GetDebugInfo(self):
        return SDP.Counter(self.pp.CFD_DEM["debug_tool_cycle"].GetInt(), 1, is_dead = 1)

    def SetBetaParameters(self):
        BaseAlgorithm.SetBetaParameters(self)

    def SetCustomBetaParameters(self, dictionary): # These are input parameters that have not yet been transferred to the interface
        var_names = [k for k in dictionary.keys()]
        var_values = [k for k in dictionary.values()]

        for name, value in zip(var_names, var_values):
            if name == 'simulation_time':
                simulation_time = value
            elif name == 'basset_force_type':
                basset_force_type = value
            elif name == 'Nq':
                Nq = value
            elif name == 'm':
                m = value
            elif name == 'number_of_quadrature_steps_in_window':
                number_of_quadrature_steps_in_window = value

        self.pp.CFD_DEM.FinalTime = simulation_time
        self.pp.CFD_DEM.fluid_already_calculated = 0
        self.pp.CFD_DEM.recovery_echo_level = 1
        self.pp.CFD_DEM.gradient_calculation_type = 5
        self.pp.CFD_DEM.pressure_grad_recovery_type = 1
        self.pp.CFD_DEM.AddEmptyValue("store_full_gradient_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("laplacian_calculation_type").SetInt(0)
        self.pp.CFD_DEM.do_search_neighbours = False
        self.pp.CFD_DEM.AddEmptyValue("material_acceleration_calculation_type").SetInt(2)
        self.pp.CFD_DEM.AddEmptyValue("faxen_force_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("vorticity_calculation_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("print_FLUID_VEL_PROJECTED_RATE_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_MATERIAL_FLUID_ACCEL_PROJECTED_option").SetBool(True)        
        self.pp.CFD_DEM.AddEmptyValue("basset_force_type").SetInt(basset_force_type)
        self.pp.CFD_DEM.AddEmptyValue("print_BASSET_FORCE_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("basset_force_integration_type").SetInt(1)
        self.pp.CFD_DEM.AddEmptyValue("n_init_basset_steps").SetInt(2)        
        self.pp.CFD_DEM.AddEmptyValue("time_steps_per_quadrature_step").SetInt(Nq)
        self.pp.CFD_DEM.AddEmptyValue("delta_time_quadrature").SetDouble( self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt() * self.pp.CFD_DEM["MaxTimeStep"].GetDouble() )
        self.pp.CFD_DEM.AddEmptyValue("quadrature_order").SetInt(2)
        self.pp.CFD_DEM.AddEmptyValue("time_window").SetDouble(0.5)
        self.pp.CFD_DEM.number_of_exponentials = m
        self.pp.CFD_DEM.number_of_quadrature_steps_in_window = number_of_quadrature_steps_in_window
        self.pp.CFD_DEM.print_steps_per_plot_step = 1
        self.pp.CFD_DEM.PostCationConcentration = False
        self.pp.CFD_DEM.AddEmptyValue("do_impose_flow_from_field_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_MATERIAL_ACCELERATION_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_VORTICITY_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_MATERIAL_ACCELERATION_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_VELOCITY_GRADIENT_option").SetBool(False)
        # Making the fluid step an exact multiple of the DEM step
        self.pp.Dt = int(self.pp.Dt / self.pp.CFD_DEM["MaxTimeStep"].GetDouble()) * self.pp.CFD_DEM["MaxTimeStep"].GetDouble()
        self.pp.viscosity_modification_type = 0.0
        self.pp.CFD_DEM.fluid_domain_volume = 0.5 ** 2 * 2 * math.pi # write down the volume you know it has

        import chandelier_parameters as ch_pp
        import chandelier as ch
        self.ch_pp = ch_pp
        self.ch = ch
        self.ch_pp.include_history_force = bool(self.pp.CFD_DEM["basset_force_type"].GetInt())

        #import quadrature as quad
        self.ch.sim = ch.AnalyticSimulator(self.ch_pp)

    def SetUpResultsDatabase(self):
        import candelier_hdf5
        self.results_database = candelier_hdf5.ResultsCandelier(self.pp, self.main_path)

    def DEMSolve(self, time = 'None'):
        self.disperse_phase_algorithm.solver.Solve()
        for node in self.spheres_model_part.Nodes:
            coor_calculated = [node.X, node.Y, node.Z]
            self.radial_error = self.results_database.CalculateError(time, coor_calculated)
            self.error_time = time

    def PerformZeroStepInitializations(self):
        # Impose initial velocity to be the terminal velocity
        self.ch.sim.CalculateNonDimensionalVars()
        terminal_velocity = self.ch.sim.NDw0 *  self.ch_pp.R *  self.ch_pp.omega

        for node in self.spheres_model_part.Nodes:
            node_to_follow_id = node.Id
            node.SetSolutionStepValue(VELOCITY_Z, terminal_velocity)

        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, self.ch_pp.u0)
            node.SetSolutionStepValue(VELOCITY_Y, self.ch_pp.v0)
            node.SetSolutionStepValue(VELOCITY_Z, 2. / 9 * 9.81 * self.ch_pp.a ** 2 / (self.ch_pp.nu * self.ch_pp.rho_f) * (self.ch_pp.rho_f - self.ch_pp.rho_p))
            node.Fix(VELOCITY_Z)
            node.SetSolutionStepValue(VELOCITY_OLD_X, self.ch_pp.u0)
            node.SetSolutionStepValue(VELOCITY_OLD_Y, self.ch_pp.v0)
            node.SetSolutionStepValue(VELOCITY_OLD_Z, 2. / 9 * 9.81 * self.ch_pp.a ** 2 / (self.ch_pp.nu * self.ch_pp.rho_f) * (self.ch_pp.rho_f - self.ch_pp.rho_p))
            node.Fix(VELOCITY_OLD_Z)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, self.ch_pp.u0)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, self.ch_pp.v0)
            node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)

    def ApplyForwardCoupling(self, alpha = 'None'):
        self.projection_module.ApplyForwardCoupling(alpha)
        for node in self.spheres_model_part.Nodes:
            x = node.X
            y = node.Y
            z = node.Z
            r = math.sqrt(x ** 2 + y ** 2)
            omega = self.ch_pp.omega
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
            new_vx = - self.ch_pp.omega * self.y
            new_vy =   self.ch_pp.omega * self.x
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
