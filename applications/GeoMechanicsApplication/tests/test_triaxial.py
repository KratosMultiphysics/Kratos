from applications.GeoMechanicsApplication.tests import test_helper
from kratos.python_scripts import KratosUnittest

class KratosGeoMechanicsTriaxialTests(KratosUnittest.TestCase):
    def test_triaxial(self):
        test_name = 'test_triaxial'
        file_path = test_helper.get_file_path(test_name)
        simulation = test_helper.run_kratos(file_path)

if __name__ == '__main__':
    KratosUnittest.main()