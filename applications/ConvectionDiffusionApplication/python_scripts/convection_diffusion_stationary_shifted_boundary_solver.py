# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver

def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionStationaryShiftedBoundarySolver(main_model_part, custom_settings)

class ConvectionDiffusionStationaryShiftedBoundarySolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    """The stationary class for convection-diffusion solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

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
            "conforming_basis" : true,
            "extension_operator_type" : "MLS",
            "mls_extension_operator_order" : 1
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        # Add heat transfer required variables
        super().AddVariables()

        # Add distance variable to represent the skin
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

    def Initialize(self):
        # Avoid zeros with positive epsilon
        tol = 1.0e-12
        for node in self.GetComputingModelPart().Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if abs(dist) < tol:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, tol)

        # Calculate the required neighbours
        nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.main_model_part)
        nodal_neighbours_process.Execute()
        elemental_neighbours_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(self.main_model_part)
        elemental_neighbours_process.Execute()

        # Create the boundary elements and MLS basis
        settings = KratosMultiphysics.Parameters("""{}""")
        settings.AddEmptyValue("model_part_name").SetString(self.main_model_part.Name)
        settings.AddEmptyValue("boundary_sub_model_part_name").SetString("shifted_boundary")
        settings.AddEmptyValue("conforming_basis").SetBool(self.settings["conforming_basis"].GetBool())
        settings.AddEmptyValue("extension_operator_type").SetString(self.settings["extension_operator_type"].GetString())
        settings.AddEmptyValue("mls_extension_operator_order").SetInt(self.settings["mls_extension_operator_order"].GetInt())
        element_type = self.settings["element_replace_settings"]["element_name"].GetString()[:-4]
        if element_type == "LaplacianShiftedBoundaryElement":
            sbm_interface_condition_name = "LaplacianShiftedBoundaryCondition"
        elif element_type == "MixedLaplacianShiftedBoundaryElement":
            sbm_interface_condition_name = "MixedLaplacianShiftedBoundaryCondition"
        else:
            raise Exception(f"Unsupported \'element_type\': {element_type}")
        settings.AddEmptyValue("sbm_interface_condition_name").SetString(sbm_interface_condition_name)
        sbm_interface_utility = KratosMultiphysics.ShiftedBoundaryMeshlessInterfaceUtility(self.model, settings)
        sbm_interface_utility.CalculateExtensionOperator()

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
        supported_elements = ["LaplacianShiftedBoundaryElement", "MixedLaplacianShiftedBoundaryElement"]
        if element_name not in supported_elements:
            raise Exception("Only \'LaplacianShiftedBoundaryElement\' and \'MixedLaplacianShiftedBoundaryElement\' are supported so far.")

        # Set registering element name
        name_string = "{0}{1}D{2}N".format(element_name, domain_size, num_nodes_elements)
        self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        # Set default conditions
        condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
        if condition_name != "":
            KratosMultiphysics.Logger.PrintWarning("ConvectionDiffusionStationaryShiftedBoundarySolver", "Ignoring provided condition \'{}\'.".format(condition_name))
        name_string = "LineCondition2D2N" if domain_size == 2 else "SurfaceCondition3D3N"
        self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        return self.settings["element_replace_settings"]