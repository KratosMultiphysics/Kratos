# Importing the Kratos Library
from numpy.lib.twodim_base import vander
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtilities

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
have_conv_diff = KratosUtilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")
if have_conv_diff:
    import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import os
# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.FluidDynamicsApplication.read_distance_from_file import DistanceImportUtility
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_two_fluids_solver import NavierStokesTwoFluidsSolver

def CreateSolver(model, custom_settings):
    return NavierStokesTwoFluidsSolverLeapFrog(model, custom_settings)

class NavierStokesTwoFluidsSolverLeapFrog(NavierStokesTwoFluidsSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {   "solver_type": "two_fluid_leap_frog",
            "time_order": 2,
            "time_scheme": "bdf2",
            "initial_first_order_steps":0,
            "maximum_spectral_radius": 0.5,
            "formulation": {
                "mass_source":true
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesTwoFluidsSolverLeapFrog, cls).GetDefaultParameters())
        return default_settings


    def __init__(self, model, custom_settings):

        super().__init__(model,custom_settings)

        # Set the time integration scheme for two fluid element.
        self.time_scheme=custom_settings["time_scheme"].GetString()

        # Crank Nicolson time integration scheme implemented as a theta scheme
        if self.time_scheme == "theta_scheme":
            self.initialize_leap_frog=True
            self.element_name= "TwoFluidNavierStokesCN"
            self.min_buffer_size = 3
            self.time_order=custom_settings["time_order"].GetInt()
            if (self.time_order == 1.0):
                self.theta = 1.0
            elif (self.time_order == 2):
                self.theta = 0.55

            else:
                raise ValueError("\'theta_scheme\' with \'time_order\' {} is not implemented. Use \'1\' for Backward Euler or \'2\' for Crank-Nicolson.".format(self.time_order))
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME_INTEGRATION_THETA, self.theta)
            self.leapfrog_coef=self.theta
        # Generalised alpha integration time scheme
        elif self.time_scheme =="generalised_alpha":
            self.element_name="TwoFluidNavierStokesAlphaMethod"
            self.min_buffer_size = 3
            self.leapfrog_coef=0.5
            self.rho_inf=custom_settings["maximum_spectral_radius"].GetInt()
            self.initialize_leap_frog=True
            self.main_model_part.ProcessInfo.SetValue(KratosCFD.SPECTRAL_RADIUS_LIMIT, self.rho_inf)

        # Bdf2 time integration time scheme. This option behaves as a TwoFLuidNavierStokesSolver.

        elif self.time_scheme == "bdf2":
            self.element_name = "TwoFluidNavierStokes"
            self.min_buffer_size = 3
        else:
            raise ValueError("\'time_scheme\' {} is not implemented. Use \'bdf2\' or \'theta_scheme\'.".format(self.time_scheme))

        #In order to avoid dragging the initial instabilities of the problem, a dissipative integration scheme such as BackWard Euler can be used for the initial time steps and then switch to CrankNicolson scheme.
        self.initial_first_order = False
        self.initial_first_order_steps = custom_settings["initial_first_order_steps"].GetInt()

        if (self.initial_first_order_steps != 0):
            if (self.time_scheme == "theta_scheme") and self.time_order == 1:
                KratosMultiphysics.Logger.PrintWarning("NavierStokesTwoFluidsSolver", "initial_first_order_steps {} for {} has not any effect since the same time scheme and order is using during simulation.".format(self.initial_first_order_steps, self.time_scheme))

            elif self.time_scheme == "bdf2":
                KratosMultiphysics.Logger.PrintWarning("NavierStokesTwoFluidsSolver", "initial_first_order_steps {} for {} has not any effect since the same time scheme and order is using during simulation.".format(self.initial_first_order_steps, self.time_scheme))

            elif self.time_scheme== "generalised_alpha":
                KratosMultiphysics.Logger.PrintWarning("NavierStokesTwoFluidsSolver", "initial_first_order_steps {} for {} has not any effect it is goint to set to 0.0.".format(self.initial_first_order_steps, self.time_scheme))
                self.initial_first_order_steps=0.0
            else:
                self.initial_first_order=True

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesTwoFluidsSolverLeapFrog finished.")


    def Initialize(self):
        super().Initialize()

        # Non historical variable is initilized in order to avoid memory problems
        if self.time_scheme == "generalised_alpha":
            for node in self.GetComputingModelPart().Nodes:
                acceleration= [0,0,0]
                node.SetValue(KratosMultiphysics.ACCELERATION,acceleration)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):

        #  Hybrid time integration scheme is activated
        if self.initial_first_order:
            time_step=self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP]
            if  time_step < self.initial_first_order_steps:
                self.theta = 1
            else:
                self.theta = 0.55
            self.leapfrog_coef=self.theta
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME_INTEGRATION_THETA, self.theta)

        # Inlet and outlet water discharge is calculated for current time step, first discharge and the considering the time step inlet and outlet volume is calculated
        if self.mass_source:
            outlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(),KratosMultiphysics.OUTLET)
            inlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(),KratosMultiphysics.INLET)
            current_dt = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            inlet_volume = -current_dt * inlet_discharge
            outlet_volume = current_dt * outlet_discharge

            # System water volume is calculated for current time step considering inlet and outlet discharge.
            system_volume = inlet_volume + self.initial_system_volume - outlet_volume

        if self._TimeBufferIsInitialized():
            if self.time_scheme == "bdf2":

                # Recompute the BDF2 coefficients
                (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

                # Perform the level-set convection according to the previous step velocity
                self._PerformLevelSetConvection()

            elif self.time_scheme == "theta_scheme" or self.time_scheme =="generalised_alpha":

            # Convection problem is moved into middle time step in order to perfom leapfrog
                if self.initialize_leap_frog:
                    # Calculate time step according to leapfrog coefficient
                    dt = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME]
                    transport_dt = self.leapfrog_coef * dt
                    self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, transport_dt)
                    # Perform the level-set convection according to the previous step velocity
                    self._PerformLevelSetConvection()
                    for node in self.GetComputingModelPart().Nodes:
                        phi_n_05 = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE, 0)
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1, phi_n_05)
                    #Restoring the user time_step
                    self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, dt)
                    # Deactivating leap_frog since after this level set convection problem the leap-frog is automatised.
                    self.initialize_leap_frog = False
                else:
                    self._PerformLevelSetConvection()

            else:
                raise Exception("Unsupported time scheme \'{}\'.".format(self.time_scheme))

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Level-set convection is performed.")

        # Perform distance correction to prevent ill-conditioned cuts
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
        self._SetNodalProperties()

        # Initialize the solver current step
        self._GetSolutionStrategy().InitializeSolutionStep()

        # Recompute the distance field according to the new level-set position
        if self._TimeBufferIsInitialized():
            if (self._reinitialization_type == "variational"):
                self._GetDistanceReinitializationProcess().Execute()
            elif (self._reinitialization_type == "parallel"):
                adjusting_parameter = 0.05
                layers = int(adjusting_parameter*self.main_model_part.GetCommunicator().GlobalNumberOfElements()) # this parameter is essential
                max_distance = 1.0 # use this parameter to define the redistancing range
                # if using CalculateInterfacePreservingDistances(), the initial interface is preserved
                self._GetDistanceReinitializationProcess().CalculateDistances(
                    self.main_model_part,
                    self._levelset_variable,
                    KratosMultiphysics.NODAL_AREA,
                    layers,
                    max_distance,
                    self._GetDistanceReinitializationProcess().CALCULATE_EXACT_DISTANCES_TO_PLANE) #NOT_CALCULATE_EXACT_DISTANCES_TO_PLANE)

            if (self._reinitialization_type != "none"):
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing process is finished.")
            # Prepare distance correction for next step
            self._GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()
        # Moving the level set with an linear extrapolation in order to enter to the to the navier stokes equation in the same time step.
        # TODO: The midle time step level set is saved in a nonhystorical variable in this case TEMPERATURE. It should be saved in other non-hystorical variable
        if self.initialize_leap_frog:
            aux = 1.0
        else:
            aux = 0.5
        for node in self.GetComputingModelPart().Nodes:
            distance_0=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            node.SetValue(KratosMultiphysics.TEMPERATURE,distance_0)
            distance_1=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE,1)
            new_distance=(distance_0-distance_1)*aux+distance_0
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,new_distance)

        # Accumulative water volume error ratio due to level set. Adding source term
        if self.mass_source:
            water_volume_after_transport = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(self.GetComputingModelPart())
            volume_error = (water_volume_after_transport - system_volume) / system_volume
            self.initial_system_volume = system_volume
        else:
            volume_error=0

        self.main_model_part.ProcessInfo.SetValue(KratosCFD.VOLUME_ERROR, volume_error)


        if self.mass_source:
            file_name="testing.txt"
            self._ExportingMassConservationData(file_name,volume_error,water_volume_after_transport,inlet_volume,outlet_volume,system_volume)

    def FinalizeSolutionStep(self):

        if self._TimeBufferIsInitialized():

            # Finalize the solver current step
            self._GetSolutionStrategy().FinalizeSolutionStep()
            # After NS equation is solverd the level set is moved again into midle time step.
            for node in self.GetComputingModelPart().Nodes:
                distance_05=node.GetValue(KratosMultiphysics.TEMPERATURE)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance_05)

            # Limit the obtained acceleration for the next step
            # This limitation should be called on the second solution step onwards (e.g. STEP=3 for BDF2)
            # We intentionally avoid correcting the acceleration in the first resolution step as this might cause problems with zero initial conditions
            if self._apply_acceleration_limitation and self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.min_buffer_size:
                self._GetAccelerationLimitationUtility().Execute()

            # According to generalised alpha method a previus time step acceleration needs to be saved.
            if self.time_scheme == "generalised_alpha":
                # Previous time step acceleration is saved an non historical variable for the next step
                alpha_m=0.5*((3-self.rho_inf)/(1+self.rho_inf))
                alpha_f=1/(1+self.rho_inf)
                gamma= 0.5+alpha_m-alpha_f
                dt=self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                for node in self.GetComputingModelPart().Nodes:
                    vn=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
                    v=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                    acceleration_n=node.GetValue(KratosMultiphysics.ACCELERATION)
                    acceleration_alpha_method=(v-vn)/(gamma*dt)+((gamma-1)/gamma)*acceleration_n
                    node.SetValue(KratosMultiphysics.ACCELERATION,acceleration_alpha_method)



    def _CreateScheme(self):
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
                domain_size,
                domain_size + 1)

            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            time_scheme = self.settings["time_scheme"].GetString()
            if (time_scheme == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                if not (time_scheme == "generalised_alpha" or time_scheme == "theta_scheme"):
                    err_msg = "Requested elemental time scheme \"" + time_scheme + "\" is not available.\n"
                    err_msg += "Available options are: \'bdf2\', \'generalised_alpha\' and \'theta_scheme\'."
                    raise Exception(err_msg)

            return scheme

    def _ExportingMassConservationData(self,filename,new_error,water_volume_after_correction,inlet,outlet,system_volume):
        if os.path.exists(filename):
            append_write = 'a' # append if already exists
        else:
            append_write = 'w' # make a new file if not
        Time=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        highscore = open(filename,append_write)
        highscore.write(str(Time)+"\t" +str(new_error) +"\t" +str(water_volume_after_correction)+"\t" +str(inlet) +'\t\t'+ str(outlet) +'\t\t'+str(new_error) +'\t\t'+str(system_volume) +'\t\t'+'\n')
        highscore.close()