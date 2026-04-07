# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ViscosityModulatorApplication as ViscosityModulatorApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Import base class file
from KratosMultiphysics.ViscosityModulatorApplication import viscosity_modulator_solver

def CreateSolver(main_model_part, custom_settings):
    return ViscosityModulatorStationarySolver(main_model_part, custom_settings)

class ViscosityModulatorStationarySolver(viscosity_modulator_solver.ViscosityModulatorSolver):
    """The stationary class for viscosity modulator solvers.

    Public member variables:
    stationary_settings -- settings for the implicit static solvers.

    See viscosity_modulator_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):

        # Construct the base solver and validate the remaining settings in the base class
        super().__init__(main_model_part, custom_settings)

        # Overwrite the base solver minimum buffer size
        buffer_2_elems = ["EulerianConvDiff","AxisymmetricEulerianConvectionDiffusion2D3N","AxisymmetricEulerianConvectionDiffusion2D4N"] #TODO: Find a better solution
        if self.settings["element_replace_settings"]["element_name"].GetString() in buffer_2_elems:
            self.min_buffer_size = 2
        else:
            self.min_buffer_size = 1

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    #### Private functions ####
    def _CreateScheme(self):
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = 1.0
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 1.0
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STATIONARY] = True

        shock_capturing_settings = self.settings["shock_capturing_settings"]
        sc_intensity = shock_capturing_settings["shock_capturing_intensity"].GetDouble()
        anisotropic_diffusion = shock_capturing_settings["use_anisotropic_diffusion"].GetBool()
        self.GetComputingModelPart().ProcessInfo.SetValue(ViscosityModulatorApplication.SHOCK_CAPTURING_INTENSITY, sc_intensity)
        self.GetComputingModelPart().ProcessInfo.SetValue(ViscosityModulatorApplication.USE_ANISOTROPIC_DISC_CAPTURING, anisotropic_diffusion)

        # As the (no) time integration is managed by the element, we set a "fake" scheme to perform the solution update
        if not self.main_model_part.IsDistributed():
            viscosity_modulator_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        else:
            viscosity_modulator_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticScheme()

        return viscosity_modulator_scheme
