from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_mechanics_solver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return EigenSolver(main_model_part, custom_settings)


class EigenSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics eigen solver.

    This class creates the mechanical solvers for eigenvalue analysis.
    It currently supports the Feast solver.

    Member variables from base class:
    settings -- Kratos parameters containing general solver settings.
    main_model_part -- the model part used to construct the solver.

    Additional member variables:
    eigensolver_settings -- settings for the eigenvalue solvers.
    """
    def __init__(self, main_model_part, custom_settings):
        settings = custom_settings.Clone()
        if settings.Has("eigensolver_settings"):
            self.eigensolver_settings = settings["eigensolver_settings"].Clone()
            settings.RemoveValue("eigensolver_settings")
        else:
            self.eigensolver_settings = KratosMultiphysics.Parameters("{}")

        # Construct the base solver.
        super().__init__(main_model_part, settings)

        # Set defaults and validate custom settings.
        default_eigensolver_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "FEAST",
            "print_feast_output": true,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "search_dimension": 10,
            "linear_solver_settings": {
                "solver_type": "skyline_lu"
            }
        }
        """)
        self.eigensolver_settings.ValidateAndAssignDefaults(default_eigensolver_settings)
        print("::[EigenSolver]:: Construction finished")

    def Initialize(self):
        print("::[EigenSolver]:: Initializing ...")
        # The solver is created here.
        super().Initialize()
        print("::[EigenSolver]:: Finished initialization.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z,self.main_model_part)
        if self.settings["pressure_dofs"].GetBool():    
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE_REACTION,self.main_model_part)
        print("::[Structural EigenSolver]:: DOF's ADDED")

    #### Private functions ####

    def _create_solution_scheme(self):
        """Create the scheme for the eigenvalue problem.

        The scheme determines the left- and right-hand side matrices in the
        generalized eigenvalue problem. 
        """
        if self.settings["solution_type"].GetString() == "Dynamic":
            solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else:
            raise Exception("Unsupported solution_type: " + self.settings["solution_type"])
        return solution_scheme

    def _create_linear_solver(self):
        """Create the eigensolver.
        
        This overrides the base class method and replaces the usual linear solver
        with an eigenvalue problem solver.
        """
        if self.eigensolver_settings["solver_type"].GetString() == "FEAST":
            feast_system_solver_settings = self.eigensolver_settings["linear_solver_settings"]
            if feast_system_solver_settings["solver_type"].GetString() == "skyline_lu":
                # default built-in feast system solver
                linear_solver = ExternalSolversApplication.FEASTSolver(self.eigensolver_settings)
            elif feast_system_solver_settings["solver_type"].GetString() == "pastix":
                feast_system_solver = ExternalSolversApplication.PastixComplexSolver(feast_system_solver_settings)
                linear_solver = ExternalSolversApplication.FEASTSolver(self.eigensolver_settings, feast_system_solver)
            else:
                raise Exception("Unsupported FEAST system solver_type: " + feast_system_solver_settings["solver_type"].GetString())
        else:
            raise Exception("Unsupported eigensolver solver_type: " + self.eigensolver_settings["solver_type"].GetString())
        return linear_solver

    def _create_mechanical_solver(self):
        eigen_scheme = self.get_solution_scheme() # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self.get_builder_and_solver() # The eigensolver is created here.
        computing_model_part = self.GetComputingModelPart()

        return StructuralMechanicsApplication.EigensolverStrategy(
            computing_model_part,
            eigen_scheme,
            builder_and_solver)
