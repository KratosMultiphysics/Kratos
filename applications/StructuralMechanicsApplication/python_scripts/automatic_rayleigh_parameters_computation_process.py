# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import eigen_solver_factory

import KratosMultiphysics.StructuralMechanicsApplication as SMA

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AutomaticRayleighComputationProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class AutomaticRayleighComputationProcess(KM.Process):
    """This class is used in order to compute automatically the Rayleigh damping parameters computing in first place the eigenvalues of the system

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing the settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KM.Process.__init__(self)

        # Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "help"                           :"This class is used in order to compute automatically the Rayleigh damping parameters computing in first place the eigenvalues of the system",
            "model_part_name"                : "Structure",
            "echo_level"                     : 0,
            "write_on_properties"            : true,
            "damping_ratio_0"                : 0.0,
            "damping_ratio_1"                : -1.0,
            "eigen_values_vector"            : [0.0],
            "eigen_system_settings" : {
                "solver_type"       : "eigen_eigensystem"
            }
        }
        """)

        # Setting solver settings
        if settings.Has("eigen_system_settings"):
            if not settings["eigen_system_settings"].Has("solver_type"):
              settings["eigen_system_settings"].AddValue("solver_type", default_parameters["eigen_system_settings"]["solver_type"])
        else:
            settings.AddValue("eigen_system_settings", default_parameters["eigen_system_settings"])
        solver_type = settings["eigen_system_settings"]["solver_type"].GetString()
        eigen_system_settings = self._auxiliary_eigen_settings(solver_type)
        default_parameters["eigen_system_settings"] = eigen_system_settings["eigen_system_settings"]

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        # We define the model parts
        self.model = Model
        self.main_model_part = self.model[self.settings["model_part_name"].GetString()]

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # The general damping ratios
        damping_ratio_0 = self.settings["damping_ratio_0"].GetDouble()
        damping_ratio_1 = self.settings["damping_ratio_1"].GetDouble()

        # We get the model parts which divide the problem
        current_process_info = self.main_model_part.ProcessInfo
        existing_computation = current_process_info.Has(SMA.EIGENVALUE_VECTOR)

        # Create auxiliary parameters
        compute_damping_coefficients_settings = KM.Parameters("""
        {
            "echo_level"          : 0,
            "damping_ratio_0"     : 0.0,
            "damping_ratio_1"     : -1.0,
            "eigen_values_vector" : [0.0]
        }
        """)

        # Setting custom parameters
        compute_damping_coefficients_settings["echo_level"].SetInt(self.settings["echo_level"].GetInt())
        compute_damping_coefficients_settings["damping_ratio_0"].SetDouble(damping_ratio_0)
        compute_damping_coefficients_settings["damping_ratio_1"].SetDouble(damping_ratio_1)

        # We check if the values are previously defined
        properties = self.main_model_part.GetProperties()
        for prop in properties:
            if prop.Has(SMA.SYSTEM_DAMPING_RATIO):
                self.settings["damping_ratio_0"].SetDouble(prop.GetValue(SMA.SYSTEM_DAMPING_RATIO))
                break
        for prop in properties:
            if prop.Has(SMA.SECOND_SYSTEM_DAMPING_RATIO):
                self.settings["damping_ratio_1"].SetDouble(prop.GetValue(SMA.SECOND_SYSTEM_DAMPING_RATIO))
                break

        # We have computed already the eigen values
        current_process_info = self.main_model_part.ProcessInfo
        precomputed_eigen_values = self.settings["eigen_values_vector"].GetVector()
        if len(precomputed_eigen_values) > 1:
            compute_damping_coefficients_settings["eigen_values_vector"].SetVector(precomputed_eigen_values)
        else:
            # If not computed eigen values already
            if not existing_computation:
                KM.Logger.PrintInfo("::[MechanicalSolver]::", "EIGENVALUE_VECTOR not previously computed. Computing automatically, take care")
                eigen_linear_solver = eigen_solver_factory.ConstructSolver(self.settings["eigen_system_settings"])
                builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(eigen_linear_solver)
                eigen_scheme = SMA.EigensolverDynamicScheme()
                eigen_solver = SMA.EigensolverStrategy(self.main_model_part, eigen_scheme, builder_and_solver,
                    self.mass_matrix_diagonal_value,
                    self.stiffness_matrix_diagonal_value)
                eigen_solver.Solve()

                # Setting the variable RESET_EQUATION_IDS
                current_process_info[SMA.RESET_EQUATION_IDS] = True

            eigenvalue_vector = current_process_info.GetValue(SMA.EIGENVALUE_VECTOR)
            compute_damping_coefficients_settings["eigen_values_vector"].SetVector(eigenvalue_vector)

        # We compute the coefficients
        coefficients_vector = SMA.ComputeDampingCoefficients(compute_damping_coefficients_settings)

        # We set the values
        if self.settings["write_on_properties"].GetBool():
            for prop in self.main_model_part.Properties:
                prop.SetValue(SMA.RAYLEIGH_ALPHA, coefficients_vector[0])
                if current_process_info.Has(KM.COMPUTE_LUMPED_MASS_MATRIX):
                    if current_process_info[KM.COMPUTE_LUMPED_MASS_MATRIX]:
                        prop.SetValue(SMA.RAYLEIGH_BETA, 0.0)
                    else:
                        prop.SetValue(SMA.RAYLEIGH_BETA, coefficients_vector[1])
                else:
                    prop.SetValue(SMA.RAYLEIGH_BETA, coefficients_vector[1])
        else:
            current_process_info.SetValue(SMA.RAYLEIGH_ALPHA, coefficients_vector[0])
            if current_process_info.Has(KM.COMPUTE_LUMPED_MASS_MATRIX):
                if current_process_info[KM.COMPUTE_LUMPED_MASS_MATRIX]:
                    current_process_info.SetValue(SMA.RAYLEIGH_BETA, 0.0)
                else:
                    current_process_info.SetValue(SMA.RAYLEIGH_BETA, coefficients_vector[1])
            else:
                current_process_info.SetValue(SMA.RAYLEIGH_BETA, coefficients_vector[1])

    def _auxiliary_eigen_settings(self, solver_type):
        """ This method returns the settings for the eigenvalues computations

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if solver_type == "feast":
            eigen_system_settings = KM.Parameters("""
            {
                "eigen_system_settings" : {
                    "solver_type"           : "feast",
                    "echo_level"            : 0,
                    "tolerance"             : 1e-10,
                    "symmetric"             : true,
                    "e_min"                 : 0.0,
                    "e_max"                 : 4.0e5,
                    "number_of_eigenvalues" : 2,
                    "subspace_size"         : 15
                }
            }
            """)
            self.mass_matrix_diagonal_value = 1.0
            self.stiffness_matrix_diagonal_value = -1.0
        else:
            eigen_system_settings = KM.Parameters("""
            {
                "eigen_system_settings" : {
                    "solver_type"       : "eigen_eigensystem"
                }
            }
            """)
            eigen_system_settings["eigen_system_settings"]["solver_type"].SetString(solver_type)
            self.mass_matrix_diagonal_value = 0.0
            self.stiffness_matrix_diagonal_value = 1.0
        return eigen_system_settings
