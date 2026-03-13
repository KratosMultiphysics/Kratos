# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_solver import ConvectionDiffusionSolver

from KratosMultiphysics import eigen_solver_factory
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

def CreateSolver(main_model_part, custom_settings):
    return EigenSolver(main_model_part, custom_settings)

class EigenSolver(ConvectionDiffusionSolver):
    """The structural mechanics eigen solver.

    This class creates the mechanical solvers for eigenvalue analysis.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        if custom_settings.Has("linear_solver_settings"):
            IssueDeprecationWarning('EigenSolver', '"linear_solver_settings" was specified which is not used in the EigenSolver. Use "eigensolver_settings"!')
            custom_settings.RemoveValue("linear_solver_settings")

        # Construct the base solver.
        super(EigenSolver, self).__init__(main_model_part, custom_settings)

        # Overwrite the base solver minimum buffer size
        buffer_2_elems = ["EulerianConvDiff","AxisymmetricEulerianConvectionDiffusion2D3N","AxisymmetricEulerianConvectionDiffusion2D4N"] #TODO: Find a better solution
        if self.settings["element_replace_settings"]["element_name"].GetString() in buffer_2_elems:
            self.min_buffer_size = 2
        else:
            self.min_buffer_size = 1

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"         : "dynamic",
            "compute_modal_decomposition": false,
            "eigensolver_settings" : {
                "solver_type"           : "spectra_sym_g_eigs_shift",
                "max_iteration"         : 1000,
                "number_of_eigenvalues" : 5,
                "echo_level"            : 1
            },
            "eigensolver_diagonal_values" : { }
        }""")
        base_parameters = super().GetDefaultParameters()
        base_parameters.RemoveValue("linear_solver_settings")
        this_defaults.AddMissingParameters(base_parameters)
        return this_defaults

    #### Private functions ####

    def _CreateScheme(self):
        """Create the scheme for the eigenvalue problem.

        The scheme determines the left- and right-hand side matrices in the
        generalized eigenvalue problem.
        """
        scheme_type = self.settings["scheme_type"].GetString()
        if scheme_type == "dynamic":
            solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else: # here e.g. a stability scheme could be added
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"dynamic\""
            raise Exception(err_msg)

        return solution_scheme

    def _CreateLinearSolver(self):
        """Create the eigensolver.

        This overrides the base class method and replaces the usual linear solver
        with an eigenvalue problem solver.
        """
        return eigen_solver_factory.ConstructSolver(self.settings["eigensolver_settings"])

    def _CreateSolutionStrategy(self):
        eigen_scheme = self._GetScheme() # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self._GetBuilderAndSolver() # The eigensolver is created here.
        computing_model_part = self.GetComputingModelPart()

        solver_type = self.settings["eigensolver_settings"]["solver_type"].GetString()
        if solver_type in ["eigen_eigensystem", "spectra_sym_g_eigs_shift", "spectra_g_eigs_shift"]: # TODO evaluate what has to be used for spectra
            mass_matrix_diagonal_value = 0.0
            stiffness_matrix_diagonal_value = 1.0
        elif solver_type == "feast":
            mass_matrix_diagonal_value = 1.0
            stiffness_matrix_diagonal_value = -1.0
        else:
            diag_values = self.settings["eigensolver_diagonal_values"]
            if not diag_values.Has("mass_matrix_diagonal_value") or not diag_values.Has("stiffness_matrix_diagonal_value"):
                err_msg  = 'For the used eigensolver "{}" no defaults for '.format(solver_type)
                err_msg += '"mass_matrix_diagonal_value" and "stiffness_matrix_diagonal_value" exist, '
                err_msg += 'please specify them under "eigensolver_diagonal_values"'
                raise Exception(err_msg)

            mass_matrix_diagonal_value = diag_values["mass_matrix_diagonal_value"].GetDouble()
            stiffness_matrix_diagonal_value = diag_values["stiffness_matrix_diagonal_value"].GetDouble()

        return StructuralMechanicsApplication.EigensolverStrategy(computing_model_part,
                                                                  eigen_scheme,
                                                                  builder_and_solver,
                                                                  mass_matrix_diagonal_value,
                                                                  stiffness_matrix_diagonal_value,
                                                                  self.settings["compute_modal_decomposition"].GetBool())