from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory

# Import scipy modules
import KratosMultiphysics.scipy_conversion_tools

def CreateSolver(main_model_part, custom_settings):
    return CustomScipySolver(main_model_part, custom_settings)

class CustomScipySolver(MechanicalSolver):
    """The structural mechanics custom scipy solver.

    This class creates the mechanical solvers to provide mass and stiffness matrices as scipy matrices.
    
    Derived class must override the function SolveSolutionStep. In there the Mass and Stiffness matrices
    can be obtained as scipy matrices.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super(CustomScipySolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[CustomScipySolver]:: ", "Construction finished")
        self.space = KratosMultiphysics.UblasSparseSpace()

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"         : "dynamic"
        }""")
        this_defaults.AddMissingParameters(super(CustomScipySolver, cls).GetDefaultSettings())
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
        ''' Linear solver will not be used. But eventually the solution stratgey calls the solver's clear function. 
        To avoid crashing linear solver is provided here'''
        return linear_solver_factory.CreateFastestAvailableDirectLinearSolver()

    def _create_mechanical_solution_strategy(self):
        if self.settings["block_builder"].GetBool() == True:
            warn_msg = "In case an eigenvalue problem is computed an elimantion builder shall be used to ensure boundary conditions are applied correctly!"
            KratosMultiphysics.Logger.PrintWarning("CustomScipySolver", warn_msg)

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
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.BUILD_LEVEL,1) #Mass Matrix
        scheme = self.get_mechanical_solution_strategy().GetScheme() 

        aux = self.get_mechanical_solution_strategy().GetSystemMatrix() 
        self.space.SetToZeroMatrix(aux)
       
        # Create dummy vectors
        b = self.space.CreateEmptyVectorPointer()
        self.space.ResizeVector( b, self.space.Size1(aux) )
        self.space.SetToZeroVector(b)

        xD = self.space.CreateEmptyVectorPointer()
        self.space.ResizeVector( xD, self.space.Size1(aux) )
        self.space.SetToZeroVector(xD)

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
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.BUILD_LEVEL,2) #Stiffness Matrix
        scheme = self.get_mechanical_solution_strategy().GetScheme() 
        
        aux = self.get_mechanical_solution_strategy().GetSystemMatrix() 
        self.space.SetToZeroMatrix(aux)

        # Create dummy vectors
        b = self.space.CreateEmptyVectorPointer()
        self.space.ResizeVector( b, self.space.Size1(aux) )
        self.space.SetToZeroVector(b)

        xD = self.space.CreateEmptyVectorPointer()
        self.space.ResizeVector( xD, self.space.Size1(aux) )
        self.space.SetToZeroVector(xD)

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
        
    def SolveSolutionStep(self):
        ## Obtain scipy matrices
        # M = self._MassMatrixComputation()
        # K = self._StiffnessMatrixComputation()

        ## Obtain Dofs
        # dofs = self.get_builder_and_solver().GetDofSet()

        return True #converged
