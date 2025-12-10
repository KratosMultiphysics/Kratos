# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import eigen_solver_factory

import KratosMultiphysics.IgaApplication as IGA
import math


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NitscheStabilizationProcess(model, settings["Parameters"])

class NitscheStabilizationProcess(KratosMultiphysics.Process):
    """This class is used in order to compute automatically the Nitsche stabilization factor.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the container of the different model parts.
    params -- Kratos parameters containing the settings.
    """
    def __init__(self, model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        params -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_condition_name" : "",
            "eigen_system_settings" : {
                "solver_type"           : "feast"
            },
            "number_of_conditions" : 1
        }""")

        # Setting solver settings
        solver_type = params["eigen_system_settings"]["solver_type"].GetString()
        eigen_system_settings = self._auxiliar_eigen_settings(solver_type)
        default_parameters["eigen_system_settings"] = eigen_system_settings["eigen_system_settings"]

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.model_part_condition = self.model[self.params["model_part_condition_name"].GetString()]

        # Call the Nitsche stabilization model Part process
        KratosMultiphysics.Process.__init__(self)
        self.process = IGA.NitscheStabilizationModelPartProcess(self.model_part_condition)
        self.process.ExecuteInitialize()

        self.model_part = self.model.GetModelPart("IgaModelPart").GetSubModelPart("Nitsche_Stabilization_" + self.params["model_part_condition_name"].GetString()[13:])

        # Define the eigenvalue size for FEAST solver
        eigenvalue_nitsche_stabilization_size = self.model_part.ProcessInfo.GetValue(IGA.EIGENVALUE_NITSCHE_STABILIZATION_SIZE)

        if solver_type == "feast":
            self.params["eigen_system_settings"]["subspace_size"].SetInt(eigenvalue_nitsche_stabilization_size)
        elif solver_type == "eigen_eigensystem":
            self.params["eigen_system_settings"]["number_of_eigenvalues"].SetInt(15)
        elif solver_type == "spectra_sym_g_eigs_shift":
            self.params["eigen_system_settings"]["number_of_eigenvalues"].SetInt(math.ceil(eigenvalue_nitsche_stabilization_size-1))

    def ExecuteInitializeSolutionStep(self):
        # Get the model parts which divide the problem
        current_process_info = self.model_part.ProcessInfo

        # Compute the eigen values
        eigen_linear_solver = eigen_solver_factory.ConstructSolver(self.params["eigen_system_settings"])
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(eigen_linear_solver)
        eigen_scheme = IGA.EigensolverNitscheStabilizationScheme()
        eigen_solver = IGA.EigensolverNitscheStabilizationStrategy(self.model_part, eigen_scheme, builder_and_solver)
        eigen_solver.Solve()

        # Compute the Nitsche stabilization factor
        eigenvalue_nitsche_stabilization_vector = current_process_info.GetValue(IGA.EIGENVALUE_NITSCHE_STABILIZATION_VECTOR)
        nitsche_stabilization_factor= eigenvalue_nitsche_stabilization_vector[eigenvalue_nitsche_stabilization_vector.Size()-1]*4*self.params["number_of_conditions"].GetInt()

        eigenvalue_nitsche_stabilization_rotation_vector = current_process_info.GetValue(IGA.EIGENVALUE_NITSCHE_STABILIZATION_ROTATION_VECTOR)
        nitsche_stabilization_rotation_factor= eigenvalue_nitsche_stabilization_rotation_vector[eigenvalue_nitsche_stabilization_rotation_vector.Size()-1]*4*self.params["number_of_conditions"].GetInt()

        # Set the Nitsche stabilization factor
        for prop in self.model_part_condition.Properties:
            prop.SetValue(IGA.NITSCHE_STABILIZATION_FACTOR, nitsche_stabilization_factor)
            prop.SetValue(IGA.NITSCHE_STABILIZATION_ROTATION_FACTOR, nitsche_stabilization_rotation_factor)

        # Reset BUILD_LEVEL to calculate the continuity enforcement matrix in coupling Nitsche condition
        self.model_part_condition.ProcessInfo.SetValue(IGA.BUILD_LEVEL,0)

    def _auxiliar_eigen_settings(self, solver_type):
        """ This method returns the settings for the eigenvalues computations

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if solver_type == "feast":
            eigen_system_settings = KratosMultiphysics.Parameters("""
            {
                "eigen_system_settings" : {
                    "solver_type"           : "feast",
                    "echo_level"            : 0,
                    "tolerance"             : 1e-10,
                    "symmetric"             : true,
                    "e_min"                 : 0.0,
                    "e_max"                 : 1.0e20,
                    "number_of_eigenvalues" : 1,
                    "subspace_size"         : 1
                }
            }
            """)
        elif solver_type == "eigen_eigensystem":
            eigen_system_settings = KratosMultiphysics.Parameters("""
            {
                "eigen_system_settings" : {
                    "solver_type"       : "eigen_eigensystem",
                    "max_iteration"         : 1000,
                    "tolerance"             : 1e-9,
                    "number_of_eigenvalues" : 2,
                    "echo_level"            : 4
                }
            }
            """)
            # eigen_system_settings["eigen_system_settings"]["solver_type"].SetString(solver_type)
        elif solver_type == "spectra_sym_g_eigs_shift":
            eigen_system_settings = KratosMultiphysics.Parameters("""
            {
                "eigen_system_settings" : {
                    "solver_type"       : "spectra_sym_g_eigs_shift",
                    "number_of_eigenvalues": 3,
                    "max_iteration": 1000,
                    "echo_level": 4
                }
            }
            """)
        return eigen_system_settings