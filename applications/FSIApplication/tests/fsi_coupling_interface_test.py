import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.FSIApplication import fsi_coupling_interface
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory

try:
    import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
    structural_available = True
except ImportError:
    structural_available = False

class FSICouplingInterfaceTest(UnitTest.TestCase):
    def test_fsi_coupling_interface(self):
        # StructuralMechanicsApplication is required since we test the Mok benchmark .mdpa
        # Note that this benchmark case requires to import the TotalLagrangianElement
        if not structural_available:
            self.skipTest("Missing required application: StructuralMechanicsApplication")
        self.print_output = False
        self.check_tolerance = 1.0e-10
        self.print_reference_values = False
        self.work_folder = "FSIProblemEmulatorTest"
        self.reference_file = "reference_embedded_symbolic_navier_stokes"
        self.setUp()
        self.runTest()
        self.tearDown()
        self.checkResults()

    def setUp(self):
        # Create a complete model part to retrieve an interface from it
        self.model = KratosMultiphysics.Model()
        model_part = self.model.CreateModelPart("MainModelPart")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Import model part from mdpa file.
        work_folder = "FSIProblemEmulatorTest"
        input_filename = 'test_FSI_emulator_Structural'
        with UnitTest.WorkFolderScope(work_folder,__file__):
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(model_part)

        # Create the convergence accelerator
        conv_acc_factory = KratosMultiphysics.Parameters("""
        {
            "solver_type": "constant_relaxation",
            "w": 0.5
        }""")
        self.convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(conv_acc_factory)

        # Create the FSI coupling interface
        fsi_coupling_int_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name": "FSICouplingInterfaceStructure",
            "parent_model_part_name": "MainModelPart.StructureInterface2D_Solid_interface",
            "input_variable_list": ["POSITIVE_FACE_PRESSURE"],
            "output_variable_list": ["DISPLACEMENT"]
        }""")
        self.fsi_coupling_interface = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            fsi_coupling_int_settings,
            self.convergence_accelerator)

    def runTest(self):
        # Set a DISPLACEMENT to be updated in the father model part
        for node in self.fsi_coupling_interface.GetFatherModelPart().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [0.01 * node.Id, 0.0, 0.0])

        # Perform the update
        self.fsi_coupling_interface.ComputeResidualVector()
        self.fsi_coupling_interface.Update()

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('test_FSI_emulator_Structural.time')

    def checkResults(self):
        # Check the previously performed update
        for node in self.fsi_coupling_interface.GetInterfaceModelPart().Nodes:
            disp_relax = node.GetSolutionStepValue(KratosMultiphysics.RELAXED_DISPLACEMENT)
            expected_disp_relax = [0.005 * node.Id, 0.0, 0.0]
            for i in range(3):
                self.assertAlmostEqual(disp_relax[i], expected_disp_relax[i])

if __name__ == '__main__':
    test = FSICouplingInterfaceTest()
    test.test_fsi_coupling_interface()