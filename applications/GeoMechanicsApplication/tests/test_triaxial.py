import os

from applications.GeoMechanicsApplication.python_scripts.gid_output_file_reader import GiDOutputFileReader
from applications.GeoMechanicsApplication.tests import test_helper
from kratos.python_scripts import KratosUnittest

class KratosGeoMechanicsTriaxialTests(KratosUnittest.TestCase):
    """
    This class contains a (regression) test for the outcomes of the triaxial experiment.
    """
    def test_triaxial(self):
        """
        The output for the test was generated for the Mohr Coulomb model using the following material properties:
        "YOUNG_MODULUS": 20000.0,
        "POISSON_RATIO": 0.25,
        "COHESION": 2.0,
        "FRICTION_ANGLE": 25.0,
        "DILATANCY_ANGLE": 2.0
        """
        test_name = 'test_triaxial'
        file_path = test_helper.get_file_path(test_name)
        simulation = test_helper.run_kratos(file_path)

        # read the output files from the simulation for comparison
        reader = GiDOutputFileReader()
        result = reader.read_output_from(os.path.join(file_path, 'output.post.res'))

        displacement = reader.nodal_values_at_time("DISPLACEMENT", 1, result)
        self.assertEqual(displacement,[[0.0, -0.2, 0.0], [0.0527793, -0.2, 0.0], [0.0, -0.100033, 0.0], [0.0524016, -0.0996909, 0.0], [0.0, 0.0, 0.0], [0.105199, -0.2, 0.0], [0.105115, -0.10005, 0.0], [0.05244, 0.0, 0.0], [0.104629, 0.0, 0.0]], "The displacement in one of the nodes is not correct.")

        stress = reader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", 1, result, [1], [0])
        self.assertEqual(stress, [[[-99.9801, -252.62, -99.9799, 0.203538, 0.0, 0.0]]], "The Cauchy stress calculation is not OK.")

        strain = reader.element_integration_point_values_at_time("ENGINEERING_STRAIN_TENSOR", 1, result, [1], [0])
        self.assertEqual(strain, [[[0.104861, -0.199728, 0.104945, 0.00044325, 0.0, 0.0]]], "The engineering strain calculation is not OK.")

if __name__ == '__main__':
    KratosUnittest.main()