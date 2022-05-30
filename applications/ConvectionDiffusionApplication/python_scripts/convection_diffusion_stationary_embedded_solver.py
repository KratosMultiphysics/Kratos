# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver

def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionStationaryEmbeddedSolver(main_model_part, custom_settings)

class ConvectionDiffusionStationaryEmbeddedSolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    """The embedded class for a stationary convection-diffusion solvers.

    Public member variables:

    See convection_diffusion_solver.py for more information.
    """

    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver and validate the remaining settings in the base class
        super().__init__(main_model_part, custom_settings)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r"""{
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        # Add heat transfer required variables
        super().AddVariables()

        # Add distance variable to represent the embedded skin
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)


    def Initialize(self):

        # Avoid level set zeros using positive epsilon
        tol = 1.0e-12
        for node in self.GetComputingModelPart().Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if abs(dist) < tol:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, tol)

        # Deactivate elements in negative distance region
        for element in self.GetComputingModelPart().Elements:
            n_pos = 0
            for node in element.GetGeometry():
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0:
                    n_pos += 1
            if n_pos == 0:
                element.Set(KratosMultiphysics.ACTIVE, False)
                for node in element.GetGeometry():
                    node.Set(KratosMultiphysics.ACTIVE, False)

        # Initialize base solver strategy
        super().Initialize()

    def _get_element_condition_replace_settings(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in (2,3):
            raise Exception("DOMAIN_SIZE not set")

        # Get element data
        num_nodes_elements = 0
        if (len(self.main_model_part.Elements) > 0):
            for elem in self.main_model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break
        num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
        element_name = self.settings["element_replace_settings"]["element_name"].GetString()

        # Element checks
        if num_nodes_elements not in (3,4):
            raise Exception("Only simplex elements are supported so far.")
        if element_name != "EmbeddedLaplacianElement":
            raise Exception("Only \'EmbeddedLaplacianElement\' is supported so far.")

        # Set registering element name
        name_string = "{0}{1}D{2}N".format(element_name, domain_size, num_nodes_elements)
        self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        # Conditions
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
