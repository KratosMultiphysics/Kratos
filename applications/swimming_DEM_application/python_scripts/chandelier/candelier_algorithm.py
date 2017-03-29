import swimming_DEM_procedures as SDP
import swimming_DEM_algorithm
BaseAlgorithm = swimming_DEM_algorithm.Algorithm

class Algorithm(BaseAlgorithm):
    def __init__(self, pp):
        BaseAlgorithm.__init__(self, pp)

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead = True)

    def GetEmbeddedCounter(self):
        return SDP.Counter(1, 3, self.pp.CFD_DEM.embedded_option)  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 4, 0)

    def GetDebugInfo(self):
        return SDP.Counter(self.pp.CFD_DEM.debug_tool_cycle, 1, is_dead = 1)

    def SetBetaParamters(self):
        pass

    def SetCustomBetaParamters(self, simulation_time, basset_force_type, Nq, m, number_of_quadrature_steps_in_window): # These are input parameters that have not yet been transferred to the interface
        self.pp.CFD_DEM.FinalTime = simulation_time
        self.pp.CFD_DEM.fluid_already_calculated = 0
        self.pp.CFD_DEM.recovery_echo_level = 1
        self.pp.CFD_DEM.gradient_calculation_type = 5
        self.pp.CFD_DEM.pressure_grad_recovery_type = 1
        self.pp.CFD_DEM.store_full_gradient = 1
        self.pp.CFD_DEM.laplacian_calculation_type = 0
        self.pp.CFD_DEM.do_search_neighbours = False
        self.pp.CFD_DEM.material_acceleration_calculation_type = 2
        self.pp.CFD_DEM.faxen_force_type = 0
        self.pp.CFD_DEM.vorticity_calculation_type = 0
        self.pp.CFD_DEM.print_FLUID_VEL_PROJECTED_RATE_option = 0
        self.pp.CFD_DEM.print_MATERIAL_FLUID_ACCEL_PROJECTED_option = True
        self.pp.CFD_DEM.basset_force_type = basset_force_type
        self.pp.CFD_DEM.print_BASSET_FORCE_option = 1
        self.pp.CFD_DEM.basset_force_integration_type = 1
        self.pp.CFD_DEM.n_init_basset_steps = 2
        self.pp.CFD_DEM.time_steps_per_quadrature_step = Nq
        self.pp.CFD_DEM.delta_time_quadrature = self.pp.CFD_DEM.time_steps_per_quadrature_step * self.pp.CFD_DEM.MaxTimeStep
        self.pp.CFD_DEM.quadrature_order = 2
        self.pp.CFD_DEM.time_window = 0.5
        self.pp.CFD_DEM.number_of_exponentials = m
        self.pp.CFD_DEM.number_of_quadrature_steps_in_window = number_of_quadrature_steps_in_window
        self.pp.CFD_DEM.print_steps_per_plot_step = 1
        self.pp.CFD_DEM.PostCationConcentration = False
        self.pp.CFD_DEM.do_impose_flow_from_field = False
        self.pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
        self.pp.CFD_DEM.print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option = False
        self.pp.CFD_DEM.print_VORTICITY_option = 0
        self.pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = False
        self.pp.CFD_DEM.print_VELOCITY_GRADIENT_option = 0
        # Making the fluid step an exact multiple of the DEM step
        self.pp.Dt = int(self.pp.Dt / self.pp.CFD_DEM.MaxTimeStep) * self.pp.CFD_DEM.MaxTimeStep
        self.pp.viscosity_modification_type = 0.0
        from math import pi
        self.pp.CFD_DEM.fluid_domain_volume = 0.5 ** 2 * 2 * pi # write down the volume you know it has
