
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
            convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        else:
            convection_diffusion_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticScheme()

        return convection_diffusion_scheme
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        # self.printDofsAndCPs() # Print control points and dofs
    

    def printDofsAndCPs(self) :
        import numpy as np
        import matplotlib.pyplot as plt
        import os

        free_node_x = []
        free_node_y = []
        free_node_z = []
        fixed_node_x = []
        fixed_node_y = []
        fixed_node_z = []
        dof = []

        z_ref = 1.0
        # Set Free the active ones
        if os.path.exists("txt_files/Id_active_control_points.txt"):
            with open('txt_files/Id_active_control_points.txt', 'r') as file:
                lines = file.readlines()
            for line in lines:
                numbers = line.split()
                node = self.main_model_part.GetNode(int(numbers[0]))
                node.Free(KratosMultiphysics.TEMPERATURE)
                node.Set(KratosMultiphysics.VISITED, False)


                if (node.Z > 0.2): continue 
                free_node_x.append(node.X)
                free_node_y.append(node.Y)
                free_node_z.append(node.Z)
                dof.append(numbers[1])
        
        dof2 = []
        free_node_x2 = []
        free_node_y2= []
        free_node_z2= []
        # Set Free the active ones
        if os.path.exists("txt_files/Id_active_control_points_condition.txt"):
            with open('txt_files/Id_active_control_points_condition.txt', 'r') as file:
                lines = file.readlines()
            for line in lines:
                numbers = line.split()
                node = self.main_model_part.GetNode(int(numbers[0]))

                if (node.Z > 0.2): continue 
                free_node_x2.append(node.X)
                free_node_y2.append(node.Y)
                free_node_z2.append(node.Z)
                dof2.append(numbers[1])

        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Cerchi grandi verdi
        ax.scatter(free_node_x2, free_node_y2, free_node_z2, marker='o', color='green', s=200, label='Free Nodes')

        # Croci rosse
        ax.scatter(free_node_x, free_node_y, free_node_z, marker='x', color='red', s=50, label='Free Nodes')


        plt.show()
