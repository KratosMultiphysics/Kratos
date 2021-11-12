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

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 1

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters(r"""{
            "mls_extension_operator_order" : 1,
            "use_mls_constraints" : true,
            "use_distance_modification" : false
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        # Add heat transfer required variables
        super().AddVariables()

        # Add distance variable to represent the embedded skin
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

    def Initialize(self):

        if self.settings["use_mls_constraints"].GetBool():
            if self.settings["use_distance_modification"].GetBool():
                KratosMultiphysics.Logger.PrintWarning("ConvectionDiffusionStationaryEmbeddedSolver", "Giving precedence to MLS constraints over distance modification")

            # Calculate the required neighbours
            nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.main_model_part)
            nodal_neighbours_process.Execute()
            avg_num_elements = 10
            dimensions = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
            elemental_neighbours_process = KratosMultiphysics.FindElementalNeighboursProcess(self.main_model_part, dimensions, avg_num_elements)
            elemental_neighbours_process.Execute()

            # Create the MLS basis and add multi point constraints for all negative nodes of intersected elements
            # NOTE: the nodal distances will still be modified slightly to avoid levelset zeros (tol=1e-12) 
            # NOTE: elements in the negative distance region will be deactivated
            settings = KratosMultiphysics.Parameters("""{}""")
            settings.AddEmptyValue("model_part_name").SetString(self.main_model_part.Name)
            settings.AddEmptyValue("mls_extension_operator_order").SetInt(self.settings["mls_extension_operator_order"].GetInt())
            settings.AddEmptyValue("avoid_zero_distances").SetBool(True)
            settings.AddEmptyValue("deactivate_negative_elements").SetBool(True)
            settings.AddEmptyValue("deactivate_intersected_elements").SetBool(False)
            constraint_process = ConvectionDiffusionApplication.EmbeddedMLSConstraintProcess(self.model, settings)
            constraint_process.Execute()

        elif self.settings["use_distance_modification"].GetBool():
            # FUTURE: use DistanceModificationProcess with tol_d
            # Avoid small cuts
            tol_d = 0.001
            for node in self.GetComputingModelPart().Nodes:
                d = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if abs(d) < tol_d:
                    if d > 0.0:
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, tol_d)
                    else:
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -tol_d)

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

        else:
            # FUTURE: use DistanceModificationProcess with tol_d
            # Avoid zero distances
            tol_d = 1.0e-12
            for node in self.GetComputingModelPart().Nodes:
                d = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if abs(d) < tol_d:
                    if d > 0.0:
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, tol_d)
                    else:
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -tol_d)

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
