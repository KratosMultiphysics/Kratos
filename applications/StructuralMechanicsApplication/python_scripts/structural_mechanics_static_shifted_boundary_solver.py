# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_static_solver

def CreateSolver(main_model_part, custom_settings):
    return StaticMechanicalShiftedBoundarySolver(main_model_part, custom_settings)

class StaticMechanicalShiftedBoundarySolver(structural_mechanics_static_solver.StaticMechanicalSolver):

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
        # Add structural mechanics required variables
        super().AddVariables()

        # Add DISTANCE variable to represent the skin
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
        settings.AddEmptyValue("sbm_interface_condition_name").SetString("DisplacementShiftedBoundaryCondition")
        sbm_interface_utility = KratosMultiphysics.ShiftedBoundaryMeshlessInterfaceUtility(self.model, settings)
        sbm_interface_utility.CalculateExtensionOperator()

        # Initialize base solver strategy
        super().Initialize()