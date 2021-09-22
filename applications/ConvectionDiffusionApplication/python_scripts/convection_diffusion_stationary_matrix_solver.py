# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver
import KratosMultiphysics.scipy_conversion_tools
import numpy as np


def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionStationaryMatrixSolver(main_model_part, custom_settings)


class ConvectionDiffusionStationaryMatrixSolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):

    """Variant of the stationary convection diffusion solver that extracts:
    - the system matrix as scipy.sparse.csr_matrix
    - the system vector as np.ndarray
    """

    def __init__(self, main_model_part, custom_settings):

        # Construct the base solver and validate the remaining settings in the base class
        super(ConvectionDiffusionStationaryMatrixSolver, self).__init__(main_model_part, custom_settings)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")


    def SolveSolutionStep(self):
        """Assembles the system and stores the system matrix and vector as member variables."""
        self.K, self.p = self._SystemComputation()

        # example of a solution with the extracted matrix and vector
        # ----
        # import scipy.sparse.linalg
        # x = scipy.sparse.linalg.spsolve(self.K, self.p)

        # self._AssignSystemVector(x)
        # ----

        KratosMultiphysics.Logger.PrintInfo("::[{}]:: ".format(self.__class__.__name__), "Extracted system matrix and vector.")
        KratosMultiphysics.Logger.PrintInfo("::[{}]:: ".format(self.__class__.__name__), "No solution triggered!")

        return True

    def _SystemComputation(self):
        """Assembles the system matrix and vector and returns them as scipy.sparse.csr_matrix and np.ndarray respectively."""
        space = KratosMultiphysics.UblasSparseSpace()
        strategy = self.get_convection_diffusion_solution_strategy()
        scheme = strategy.GetScheme()

        A = strategy.GetSystemMatrix()
        space.SetToZeroMatrix(A)

        b = strategy.GetSystemVector()
        space.SetToZeroVector(b)

        # Create dummy vector
        xD = space.CreateEmptyVectorPointer()
        space.ResizeVector( xD, space.Size1(A) )
        space.SetToZeroVector(xD)

        # Build matrix
        builder_and_solver = self.get_builder_and_solver()
        builder_and_solver.Build(scheme, self.GetComputingModelPart(), A, b)
        # Apply constraints
        builder_and_solver.ApplyConstraints(scheme, self.GetComputingModelPart(), A, b)
        # Apply boundary conditions
        builder_and_solver.ApplyDirichletConditions(scheme, self.GetComputingModelPart(), A, xD, b)
        # Convert system matrix to scipy
        A_csr = KratosMultiphysics.scipy_conversion_tools.to_csr(A)
        # Convert system vector to np
        b_np = np.array(b)

        return A_csr, b_np


    def _AssignSystemVector(self, vector):
        """Assigns the values of the vector to the TEMPERATURE dofs."""
        for node in self.GetComputingModelPart().Nodes:
            dof = node.GetDof(KratosMultiphysics.TEMPERATURE)

            if not dof.IsFixed():
                value = vector[dof.EquationId]
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, value)
