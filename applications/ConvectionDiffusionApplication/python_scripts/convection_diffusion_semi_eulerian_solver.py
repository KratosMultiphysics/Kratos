from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_base_solver

def CreateSolver(model, custom_settings):
    return ConvectionDiffusionSemiEulerianSolver(model, custom_settings)

class ConvectionDiffusionSemiEulerianSolver(convection_diffusion_base_solver.ConvectionDiffusionBaseSolver):
    """The semi-eulerian class for convection-diffusion solvers.

    See convection_diffusion_base_solver.py for more information.
    """

    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the remaining settings in the base class
        super().__init__(model, custom_settings)

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "analysis_type" : "linear",
            "time_integration_method" : "implicit",
            "solver_type" : "convection_diffusion_semi_eulerian_solver",
            "element_replace_settings" : {
                "element_name" : "EulerianDiffusion",
                "condition_name" : "ThermalFace"
            },
            "reform_dofs_at_each_step" : false,
            "compute_reactions" : true,
            "pure_convection" : false,
            "bfecc_substepping": 10
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def SolveSolutionStep(self):
        # Perform convection
        bfecc_substepping = self.settings["bfecc_substepping"].GetInt()
        thermal_settings  = self.GetComputingModelPart().ProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)
        unknown_variable = thermal_settings.GetUnknownVariable()
        projection_variable = thermal_settings.GetProjectionVariable()
        velocity_variable = thermal_settings.GetVelocityVariable()

        VariableUtils().CopyScalarVar(unknown_variable, projection_variable, self.GetComputingModelPart())
        self.bfecc_utility.CopyScalarVarToPreviousTimeStep(self.GetComputingModelPart(), projection_variable)
        self.bfecc_utility.BFECCconvect(self.GetComputingModelPart(), projection_variable, velocity_variable, bfecc_substepping)

        # Solve diffusion
        is_converged = True
        pure_convection = self.settings["pure_convection"].GetBool()
        if not pure_convection:
            is_converged = self.get_convection_diffusion_solution_strategy().SolveSolutionStep()

        return is_converged

    #### Private functions ####
    def _get_element_condition_replace_settings(self):
        # Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in (2,3):
            raise Exception("DOMAIN_SIZE not set")

        # Set pure diffusion elements
        if self.settings["element_replace_settings"].Has("element_name"):
            element_name = self.settings["element_replace_settings"]["element_name"].GetString()
            warning_msg = "\'element_replace_settings\' has \'element_name\': {0}. Setting required \'EulerianDiffusion\'".format(element_name)
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, warning_msg)
        self.settings["element_replace_settings"]["element_name"].SetString("EulerianDiffusion")

        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            for elem in self.main_model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break
        num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)

        if domain_size == 2:
            if num_nodes_elements == 3:
                self.settings["element_replace_settings"]["element_name"].SetString("EulerianDiffusion2D3N")
            else:
                # TODO: Check that the EulerianDiffusion2D4N works!
                # self.settings["element_replace_settings"]["element_name"].SetString("EulerianDiffusion2D4N")
                err_msg = "EulerianDiffusion2D4N not exported yet."
                raise Exception(err_msg)
        else:
            if num_nodes_elements == 4:
                self.settings["element_replace_settings"]["element_name"].SetString("EulerianConvDiff3D4N")
            else:
                # TODO: Check that the EulerianDiffusion3D8N works!
                # self.settings["element_replace_settings"]["element_name"].SetString("EulerianDiffusion3D8N")
                err_msg = "EulerianDiffusion3D8N not exported yet."
                raise Exception(err_msg)

        # TODO: Create an auxiliary method to get this from the base
        # Set thermal conditions
        num_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().SumAll(len(self.main_model_part.Conditions))

        if num_conditions > 0:
            num_nodes_conditions = 0
            if (len(self.main_model_part.Conditions) > 0):
                for cond in self.main_model_part.Conditions:
                    num_nodes_conditions = len(cond.GetNodes())
                    break

            num_nodes_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)

            condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
            if condition_name in ("FluxCondition","ThermalFace","Condition"):
                name_string = "{0}{1}D{2}N".format(condition_name,domain_size, num_nodes_conditions)
                self.settings["element_replace_settings"]["condition_name"].SetString(name_string)
        else:
            self.settings["element_replace_settings"]["condition_name"].SetString("")

        return self.settings["element_replace_settings"]

    def _create_solution_scheme(self):
        # Create a "fake" time scheme to perform the solution update
        convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return convection_diffusion_scheme

    def _create_convection_diffusion_solution_strategy(self):
        # Create the solution strategy for the diffusion problem
        diffusion_solution_strategy = super().convection_diffusion_solution_strategy()

        # TODO: Encapsulate these in auxiliary methods
        # Create the locator for the BFECC convector
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 2:
            self.locator = BinBasedFastPointLocator2D(self.main_model_part)
        elif domain_size == 3:
            self.locator = BinBasedFastPointLocator3D(self.main_model_part)
        else:
            err_msg = "Wrong domain size: {0}".format(domain_size)
            raise Exception(err_msg)

        # Initialize the locator search database
        self.locator.UpdateSearchDatabase()

        # Create the BFECC convection utility
        if domain_size ==2:
            self.bfecc_utility = BFECCConvection2D(self.locator)
        else:
            self.bfecc_utility = BFECCConvection3D(self.locator)

        return diffusion_solution_strategy
