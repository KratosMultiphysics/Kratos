import os
import KratosMultiphysics                as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsMohrCoulombWithTensionTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """

    def simulate_mohr_coulomb(self, test_name):
        file_path = test_helper.get_file_path(os.path.join('test_mohr_coulomb_with_tension_cutoff', test_name))
        simulation = test_helper.run_kratos(file_path)
        return test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_TENSOR)

    def test_dirichlet_corner_return_zone(self):
        cauchy_stresses = self.simulate_mohr_coulomb('test_dirichlet_corner_return_zone')
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_xx = sig_integrationpoint1_element1[0,0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_yy, -1.5179192179966735)

    def test_dirichlet_elastic_zone(self):
        cauchy_stresses = self.simulate_mohr_coulomb('test_dirichlet_elastic_zone')
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_xx = sig_integrationpoint1_element1[0,0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        self.assertAlmostEqual(sig_xx, -10.0)
        self.assertAlmostEqual(sig_yy, -10.0)

    def test_dirichlet_regular_failure_zone(self):
        cauchy_stresses = self.simulate_mohr_coulomb('test_dirichlet_regular_failure_zone')
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_xx = sig_integrationpoint1_element1[0,0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        self.assertAlmostEqual(sig_xx, 3.430712712091948)
        self.assertAlmostEqual(sig_yy, -25.759721409731483)

    def test_dirichlet_tension_apex_return_zone(self):
        cauchy_stresses = self.simulate_mohr_coulomb('test_dirichlet_tension_apex_return_zone')
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_xx = sig_integrationpoint1_element1[0,0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_yy, 10.0)

    def test_dirichlet_tension_cutoff_return_zone(self):
        cauchy_stresses = self.simulate_mohr_coulomb('test_dirichlet_tension_cutoff_return_zone')
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_xx = sig_integrationpoint1_element1[0,0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_yy, 8.0)
        
if __name__ == '__main__':
    KratosUnittest.main()
