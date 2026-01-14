import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper


class KratosGeoMechanicsCaseSensitivityInUMATTests(KratosUnittest.TestCase):
    def test_case_sensitive_umats(self):
        test_name = 'test_case_sensitive_UMATS'
        project_path = test_helper.get_file_path(test_name)
        n_stages = 1
        run_multiple_stages.run_stages(project_path, n_stages)


if __name__ == '__main__':
    KratosUnittest.main()
