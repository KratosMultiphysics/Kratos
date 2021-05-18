# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver

def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionStationaryShiftedBoundarySolver(main_model_part, custom_settings)

class ConvectionDiffusionStationaryShiftedBoundarySolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    """The stationary class for convection-diffusion solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver and validate the remaining settings in the base class
        super().__init__(main_model_part, custom_settings)

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 1

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    def AddVariables(self):
        # Add heat transfer required variables
        super().AddVariables()

        # Add distance variable to represent the skin
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

    def Initialize(self):
        # Initialize base solver strategy
        super().Initialize()

        # Find the required problem subdomains
        self.__SetInterfaceFlags()

        # Calculate the nodal neighbours
        nodal_neighbours_process = KratosMultiphysics.FindNodalNeighboursProcess(self.main_model_part)
        nodal_neighbours_process.Execute()

        # Create the boundary elements and MLS basis

    def __SetInterfaceFlags(self):
        # Initialize the required flags to false
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, False, self.GetComputingModelPart().Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, False, self.GetComputingModelPart().Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, False, self.GetComputingModelPart().Elements)

        # Mark the first layer of positive distance nodes as BOUNDARY
        for element in self.GetComputingModelPart().Elements:
            geometry = element.GetGeometry()
            if self.__IsSplit(geometry):
                for node in geometry:
                    dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                    if dist > 0.0:
                        node.Set(KratosMultiphysics.BOUNDARY, True)

        # Mark the first layer of positive distance elements and their nodes as INTERFACE
        for element in self.GetComputingModelPart().Elements:
            geometry = element.GetGeometry()
            if self.__IsPositive(geometry):
                # Check if one of the nodes is BOUNDARY
                is_boundary = False
                for node in geometry:
                    if node.Is(KratosMultiphysics.BOUNDARY):
                        is_boundary = True
                        break
                # If BOUNDARY, mark the element and all the nodes as INTERFACE
                if is_boundary:
                    element.Set(KratosMultiphysics.INTERFACE, True)
                    for node in geometry:
                        node.Set(KratosMultiphysics.INTERFACE, True)

    def __IsSplit(self, geometry):
        n_neg = 0
        n_pos = 0
        for node in geometry:
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                n_neg += 1
            else:
                n_pos += 1
        return (n_neg != 0 and n_pos != 0)

    def __IsPositive(self, geometry):
        n_pos = 0
        for node in geometry:
            if not node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                n_pos += 1
        return n_pos == len(geometry)