
import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities  import ReadModelPart
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestNeighbourUtils(kratos_unittest.TestCase):
    def test_InitializeParentElementForConditions(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Structure")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        ref_map = {
            1: [1, 2],
            2: [2, 3],
            3: [2, 3],
            4: [1, 2]
        }

       # create the primal analysis execution policy wrapper
        with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
            ReadModelPart("Structure", model_part)
            KratosOA.NeighbourUtils.InitializeParentElementForConditions(model_part)
            computed_map = KratosOA.NeighbourUtils.GetConditionIdAndParentElementIdMap(model_part)
            for k, v in computed_map.items():
                self.assertEqual(ref_map[k], v)

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope("linear_strain_energy_test", __file__):
            DeleteFileIfExisting("Structure.time")

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()