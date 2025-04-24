import os
import KratosMultiphysics                as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsMohrCoulombWithTensionTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """

    def simulate_mohr_coulomb(self, test_name, dimension):
        file_path = test_helper.get_file_path(os.path.join('test_mohr_coulomb_with_tension_cutoff', test_name))
        simulation = test_helper.run_kratos(file_path)
        cauchy_stresses = test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_TENSOR)
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_xx = sig_integrationpoint1_element1[0,0]
        index = dimension - 1
        sig_yy_or_zz = sig_integrationpoint1_element1[index,index]
        return sig_xx, sig_yy_or_zz

    def test_dirichlet_corner_return_zone_2d(self):
        sig_xx, sig_yy = self.simulate_mohr_coulomb('test_dirichlet_corner_return_zone_2d', 2)
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_yy, -1.5179192179966735)

    def test_dirichlet_elastic_zone_2d(self):
        sig_xx, sig_yy = self.simulate_mohr_coulomb('test_dirichlet_elastic_zone_2d', 2)
        self.assertAlmostEqual(sig_xx, -10.0)
        self.assertAlmostEqual(sig_yy, -10.0)

    def test_dirichlet_regular_failure_zone_2d(self):
        sig_xx, sig_yy = self.simulate_mohr_coulomb('test_dirichlet_regular_failure_zone_2d', 2)
        self.assertAlmostEqual(sig_xx, 3.430712712091948)
        self.assertAlmostEqual(sig_yy, -25.759721409731483)

    def test_dirichlet_tension_apex_return_zone_2d(self):
        sig_xx, sig_yy = self.simulate_mohr_coulomb('test_dirichlet_tension_apex_return_zone_2d', 2)
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_yy, 10.0)

    def test_dirichlet_tension_cutoff_return_zone_2d(self):
        sig_xx, sig_yy = self.simulate_mohr_coulomb('test_dirichlet_tension_cutoff_return_zone_2d', 2)
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_yy, 8.0)

    def test_dirichlet_corner_return_zone_3d(self):
        sig_xx, sig_zz = self.simulate_mohr_coulomb('test_dirichlet_corner_return_zone_3d', 3)
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_zz, -1.5179192179966735)

    def test_dirichlet_elastic_zone_3d(self):
        sig_xx, sig_zz = self.simulate_mohr_coulomb('test_dirichlet_elastic_zone_3d', 3)
        self.assertAlmostEqual(sig_xx, -10.0)
        self.assertAlmostEqual(sig_zz, -10.0)

    def test_dirichlet_regular_failure_zone_3d(self):
        sig_xx, sig_zz = self.simulate_mohr_coulomb('test_dirichlet_regular_failure_zone_3d', 3)
        self.assertAlmostEqual(sig_xx, 3.430712712091948)
        self.assertAlmostEqual(sig_zz, -25.759721409731483)

    def test_dirichlet_tension_apex_return_zone_3d(self):
        sig_xx, sig_zz = self.simulate_mohr_coulomb('test_dirichlet_tension_apex_return_zone_3d', 3)
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_zz, 10.0)

    def test_dirichlet_tension_cutoff_return_zone_3d(self):
        sig_xx, sig_zz = self.simulate_mohr_coulomb('test_dirichlet_tension_cutoff_return_zone_3d', 3)
        self.assertAlmostEqual(sig_xx, 10.0)
        self.assertAlmostEqual(sig_zz, 8.0)
        
    def test_neumann_direchlet_interface_mohr_coulomb_2plus2(self):
        test_name = "test_neumann_direchlet_interface_mohr_coulomb_2plus2"
        file_path = test_helper.get_file_path(os.path.join('test_mohr_coulomb_with_tension_cutoff', test_name))
        simulation = test_helper.run_kratos(file_path)
        reaction = test_helper.get_nodal_variable(simulation, Kratos.REACTION)
        self.assertAlmostEqual(-176.1686911267587877808, reaction[0][0])
        
if __name__ == '__main__':
    KratosUnittest.main()
