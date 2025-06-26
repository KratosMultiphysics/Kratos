import os
import shutil
import importlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the Kratos applications that we need
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

import test_helper

class KratosGeoMechanicsSettlementWorkflowCppRoute(KratosUnittest.TestCase):
    """
    This test class checks the settlement workflow using the C++ route.
    """
    def get_test_dir_name(self):
        return "cpp"

    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("moving_column_with_fixed_pressure_above_phreatic_line")
        self.test_path = os.path.join(self.test_root, self.get_test_dir_name())

    def test_d_settlement_workflow(self):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        status = run_geo_settlement.run_stages(self.test_root,["ProjectParameters.json"])
        self.assertEqual(status, 0)


if __name__ == '__main__':
    KratosUnittest.main()
