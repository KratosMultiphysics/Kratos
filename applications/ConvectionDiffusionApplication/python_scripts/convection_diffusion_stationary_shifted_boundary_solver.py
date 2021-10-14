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
            "mls_extension_operator_order" : 1,
            "mls_conforming_basis" : true,
            "lagrange_multipliers_imposition" : false
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        # Add heat transfer required variables
        super().AddVariables()

        # Add distance variable to represent the skin
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # Add Lagrange multipliers required variables
        if self.settings["lagrange_multipliers_imposition"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER)

    def AddDofs(self):
        # Add heat transfer DOFs
        super().AddDofs()

        # Add Lagrange multipliers DOFs
        if self.settings["lagrange_multipliers_imposition"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER, self.main_model_part)

    def Initialize(self):
        # # Correct the level set
        # #TODO: Use the FluidDynamicsApplication process
        # #FIXME: USE THIS FOR THE ZIG-ZAG PLATE TESTS
        # tol = 1.0e-12
        # for node in self.GetComputingModelPart().Nodes:
        #     dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
        #     if abs(dist) < tol:
        #         print("Node {} dist {}".format(node.Id, dist))
        #         if dist < 0.0:
        #             node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -tol)
        #         else:
        #             node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, tol)
        # #FIXME: USE THIS FOR THE ZIG-ZAG PLATE TESTS

        # Avoid zeros with positive epsilon
        tol = 1.0e-12
        for node in self.GetComputingModelPart().Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if abs(dist) < tol:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, tol)

        # Calculate the required neighbours
        nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.main_model_part)
        nodal_neighbours_process.Execute()
        avg_num_elements = 10
        dimensions = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        elemental_neighbours_process = KratosMultiphysics.FindElementalNeighboursProcess(self.main_model_part, dimensions, avg_num_elements)
        elemental_neighbours_process.Execute()

        # Create the boundary elements and MLS basis
        settings = KratosMultiphysics.Parameters("""{}""")
        settings.AddEmptyValue("model_part_name").SetString(self.main_model_part.Name)
        settings.AddEmptyValue("boundary_sub_model_part_name").SetString("shifted_boundary")
        settings.AddEmptyValue("mls_extension_operator_order").SetInt(self.settings["mls_extension_operator_order"].GetInt())
        settings.AddEmptyValue("mls_conforming_basis").SetBool(self.settings["mls_conforming_basis"].GetBool())
        if self.settings["lagrange_multipliers_imposition"].GetBool():
            sbm_interface_condition_name = "LaplacianShiftedBoundaryLagrangeMultipliersCondition"
        else:
            sbm_interface_condition_name = "LaplacianShiftedBoundaryCondition"
        settings.AddEmptyValue("sbm_interface_condition_name").SetString(sbm_interface_condition_name)
        sbm_interface_process = ConvectionDiffusionApplication.ShiftedBoundaryMeshlessInterfaceProcess(self.model, settings)
        sbm_interface_process.Execute()

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
        if element_name != "LaplacianShiftedBoundaryElement":
            raise Exception("Only \'LaplacianShiftedBoundaryElement\' is supported so far.")

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