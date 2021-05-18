from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_solver

def CreateSolver(model, custom_settings):
    return ConvectionDiffusionSemiImplicitSolver(model, custom_settings)

class ConvectionDiffusionSemiImplicitSolver(convection_diffusion_solver.ConvectionDiffusionSolver):
    """The semi-eulerian class for convection-diffusion solvers.

    See convection_diffusion_solver.py for more information.
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
            "time_integration_method" : "semi_implicit",
            "solver_type" : "transient",
            "element_replace_settings" : {
                "element_name" : "EulerianDiffusion",
                "condition_name" : "ThermalFace"
            },
            "pure_convection" : false,
            "bfecc_substepping": 10
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def Initialize(self):
        super().Initialize()

        # Trigger the BFECC convection to create it
        self._GetBFECCConvection()

    def SolveSolutionStep(self):
        # Perform the semi-Eulerian convection
        bfecc_substepping = self.settings["bfecc_substepping"].GetInt()
        thermal_settings  = self.GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        unknown_variable = thermal_settings.GetUnknownVariable()
        projection_variable = thermal_settings.GetProjectionVariable()
        velocity_variable = thermal_settings.GetVelocityVariable()

        KratosMultiphysics.VariableUtils().CopyVariable(unknown_variable, projection_variable, self.GetComputingModelPart().Nodes)
        self._GetBFECCConvection().CopyScalarVarToPreviousTimeStep(self.GetComputingModelPart(), projection_variable)
        self._GetBFECCConvection().BFECCconvect(self.GetComputingModelPart(), projection_variable, velocity_variable, bfecc_substepping)

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
                self.settings["element_replace_settings"]["element_name"].SetString("EulerianDiffusion3D4N")
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

    def _GetBFECCConvection(self):
        if not hasattr(self, '_bfecc_convection'):
            self._bfecc_convection = self._CreateBFECCConvection()
        return self._bfecc_convection

    def _create_solution_scheme(self):
        # Create a "fake" time scheme to perform the solution update
        convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return convection_diffusion_scheme

    def _CreateBFECCConvection(self):
        # Create the locator for the BFECC convector
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 2:
            point_locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.GetComputingModelPart())
        elif domain_size == 3:
            point_locator = KratosMultiphysics.BinBasedFastPointLocator3D(self.GetComputingModelPart())
        else:
            err_msg = "Wrong domain size: {0}".format(domain_size)
            raise Exception(err_msg)

        # Initialize the locator search database
        point_locator.UpdateSearchDatabase()

        # Create the BFECC convection utility
        if domain_size ==2:
            bfecc_utility = ConvectionDiffusionApplication.BFECCConvection2D(point_locator)
        else:
            bfecc_utility = ConvectionDiffusionApplication.BFECCConvection3D(point_locator)

        return bfecc_utility

