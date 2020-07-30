from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Import scipy modules
import KratosMultiphysics.scipy_conversion_tools
from scipy.sparse.linalg import eigsh

def CreateSolver(main_model_part, custom_settings):
    return CustomScipyBaseSolver(main_model_part, custom_settings)

class CustomScipyBaseSolver(MechanicalSolver):
    """The structural mechanics custom scipy base solver.

    This class creates the mechanical solvers to provide mass and stiffness matrices as scipy matrices.

    Derived class must override the function SolveSolutionStep. In there the Mass and Stiffness matrices
    can be obtained as scipy matrices.
    The computation of the eigenvalue problem in this implementation is only an example how this solver is to be used.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super(CustomScipyBaseSolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[CustomScipyBaseSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"         : "dynamic"
        }""")
        this_defaults.AddMissingParameters(super(CustomScipyBaseSolver, cls).GetDefaultSettings())
        return this_defaults

    #### Private functions ####
    def _create_solution_scheme(self):
        """Create the scheme for the scipy solver.

        The scheme determines the mass and stiffness matrices
        """
        scheme_type = self.settings["scheme_type"].GetString()
        if scheme_type == "dynamic":
            solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else: # here e.g. a stability scheme could be added
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"dynamic\""
            raise Exception(err_msg)

        return solution_scheme

    def _create_linear_solver(self):
        ''' Linear solver will not be used. But eventually the solution strategy calls the solver's clear function.
        To avoid crashing linear solver is provided here'''
        return KratosMultiphysics.LinearSolver()

    def _create_mechanical_solution_strategy(self):
        if self.settings["block_builder"].GetBool():
            warn_msg = "In case an eigenvalue problem is computed an elimantion builder shall be used to ensure boundary conditions are applied correctly!"
            KratosMultiphysics.Logger.PrintWarning("CustomScipyBaseSolver", warn_msg)

        eigen_scheme = self.get_solution_scheme() # The scheme defines the matrices
        computing_model_part = self.GetComputingModelPart()
        builder_and_solver = self.get_builder_and_solver()
        linear_solver = self.get_linear_solver()

        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              eigen_scheme,
                                                              linear_solver,
                                                              builder_and_solver,
                                                              False,
                                                              False,
                                                              False,
                                                              False )

    def _MassMatrixComputation(self):
        space = KratosMultiphysics.UblasSparseSpace()
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.BUILD_LEVEL,1) #Mass Matrix
        scheme = self.get_mechanical_solution_strategy().GetScheme()

        aux = self.get_mechanical_solution_strategy().GetSystemMatrix()
        space.SetToZeroMatrix(aux)

        # Create dummy vectors
        b = space.CreateEmptyVectorPointer()
        space.ResizeVector( b, space.Size1(aux) )
        space.SetToZeroVector(b)

        xD = space.CreateEmptyVectorPointer()
        space.ResizeVector( xD, space.Size1(aux) )
        space.SetToZeroVector(xD)

        # Build matrix
        builder_and_solver = self.get_builder_and_solver()
        builder_and_solver.Build(scheme, self.GetComputingModelPart(), aux, b)
        # Apply Constraints
        builder_and_solver.ApplyConstraints(scheme, self.GetComputingModelPart(), aux, b)
        # Apply Boundary Conditions
        builder_and_solver.ApplyDirichletConditions(scheme, self.GetComputingModelPart(), aux, xD, b)
        # Convert Mass matrix to scipy
        M = KratosMultiphysics.scipy_conversion_tools.to_csr(aux)

        return M

    def _StiffnessMatrixComputation(self):
        space = KratosMultiphysics.UblasSparseSpace()
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.BUILD_LEVEL,2) #Stiffness Matrix
        scheme = self.get_mechanical_solution_strategy().GetScheme()

        aux = self.get_mechanical_solution_strategy().GetSystemMatrix()
        space.SetToZeroMatrix(aux)

        # Create dummy vectors
        b = space.CreateEmptyVectorPointer()
        space.ResizeVector( b, space.Size1(aux) )
        space.SetToZeroVector(b)

        xD = space.CreateEmptyVectorPointer()
        space.ResizeVector( xD, space.Size1(aux) )
        space.SetToZeroVector(xD)

        # Build matrix
        builder_and_solver = self.get_builder_and_solver()
        builder_and_solver.Build(scheme, self.GetComputingModelPart(), aux, b)
        # Apply constraints
        builder_and_solver.ApplyConstraints(scheme, self.GetComputingModelPart(), aux, b)
        # Apply boundary conditions
        builder_and_solver.ApplyDirichletConditions(scheme, self.GetComputingModelPart(), aux, xD, b)
        # Convert stiffness matrix to scipy
        K = KratosMultiphysics.scipy_conversion_tools.to_csr(aux)

        return K

    def _AssignVariables(self, eigenvalues, eigenvectors):
        num_eigenvalues = eigenvalues.size
        # Store eigenvalues in process info
        eigenvalue_vector = self.GetComputingModelPart().ProcessInfo.GetValue(StructuralMechanicsApplication.EIGENVALUE_VECTOR)
        eigenvalue_vector.Resize(num_eigenvalues)
        for i in range(num_eigenvalues):
            eigenvalue_vector[i] = eigenvalues[i]
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.EIGENVALUE_VECTOR, eigenvalue_vector)

        # Store eigenvectors in nodes
        for node in self.GetComputingModelPart().Nodes:
            node_eigenvectors = node.GetValue(StructuralMechanicsApplication.EIGENVECTOR_MATRIX)
            if self.settings["rotation_dofs"].GetBool() == True:
                dofs = [node.GetDof(KratosMultiphysics.ROTATION_X),
                        node.GetDof(KratosMultiphysics.ROTATION_Y),
                        node.GetDof(KratosMultiphysics.ROTATION_Z),
                        node.GetDof(KratosMultiphysics.DISPLACEMENT_X),
                        node.GetDof(KratosMultiphysics.DISPLACEMENT_Y),
                        node.GetDof(KratosMultiphysics.DISPLACEMENT_Z)]

                node_eigenvectors.Resize(num_eigenvalues, 6 )
            else:
                dofs = [node.GetDof(KratosMultiphysics.DISPLACEMENT_X),
                        node.GetDof(KratosMultiphysics.DISPLACEMENT_Y),
                        node.GetDof(KratosMultiphysics.DISPLACEMENT_Z)]
                node_eigenvectors.Resize(num_eigenvalues, 3 )

            # Fill the eigenvector matrix
            for i in range(num_eigenvalues):
                j = -1
                for dof in dofs:
                    j = j + 1
                    if dof.IsFixed():
                        node_eigenvectors[i,j] = 0.0
                    else:
                        node_eigenvectors[i,j] = eigenvectors[dof.EquationId,i]
            node.SetValue(StructuralMechanicsApplication.EIGENVECTOR_MATRIX, node_eigenvectors)

    def SolveSolutionStep(self):
        """This method must be overriden in derived class.

        The computation of the egenvalue problem is only an example how this solver is to be used.
        """
        ## Obtain scipy matrices
        M = self._MassMatrixComputation()
        K = self._StiffnessMatrixComputation()

        ## Compute eigenvalues and eigenvectors
        tolerance = 1e-6
        iteration = M.size*100
        vals, vecs = eigsh(K, 5, M, which='SM', tol=tolerance, maxiter = iteration)

        ## Assign results to Kratos variables
        self._AssignVariables(vals,vecs)

        return True #converged
