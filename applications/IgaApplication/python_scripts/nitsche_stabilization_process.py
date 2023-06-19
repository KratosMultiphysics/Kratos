# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import eigen_solver_factory

import KratosMultiphysics.IgaApplication as IGA


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
                    "solver_type"           : "feast",
                    "echo_level"            : 0,
                    "tolerance"             : 1e-10,
                    "symmetric"             : true,
                    "e_min"                 : 0.0,
                    "e_max"                 : 1.0e20,
                    "number_of_eigenvalues" : 1,
                    "subspace_size"         : 1
            },
            "number_of_conditions" : 1
        }""")

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
        self.params["eigen_system_settings"]["subspace_size"].SetInt(eigenvalue_nitsche_stabilization_size)
        self.params["eigen_system_settings"]["number_of_eigenvalues"].SetInt(eigenvalue_nitsche_stabilization_size)

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

        # Set the Nitsche stabilization factor
        for prop in self.model_part_condition.Properties:
            prop.SetValue(IGA.NITSCHE_STABILIZATION_FACTOR, nitsche_stabilization_factor)

        # Reset BUILD_LEVEL to calculate the continuity enforcement matrix in coupling Nitsche condition
        self.model_part_condition.ProcessInfo.SetValue(IGA.BUILD_LEVEL,0)