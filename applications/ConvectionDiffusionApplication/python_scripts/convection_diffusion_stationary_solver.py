
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_solver

def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionStationarySolver(main_model_part, custom_settings)

class ConvectionDiffusionStationarySolver(convection_diffusion_solver.ConvectionDiffusionSolver):
    """The stationary class for convection-diffusion solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):

        # Construct the base solver and validate the remaining settings in the base class
        super(ConvectionDiffusionStationarySolver, self).__init__(main_model_part, custom_settings)

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
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 0.0

        # As the (no) time integration is managed by the element, we set a "fake" scheme to perform the solution update
        if not self.main_model_part.IsDistributed():
            # convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
            scheme_settings = KratosMultiphysics.Parameters("""{}""")
            convection_diffusion_scheme = KratosMultiphysics.Future.NewScheme(
                self.GetComputingModelPart(),
                scheme_settings)
        else:
            convection_diffusion_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticScheme()

        return convection_diffusion_scheme
