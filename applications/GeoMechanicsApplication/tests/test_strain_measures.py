import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsStrainMeasureTests(KratosUnittest.TestCase):
    """
    This class contains tests which check the displacement result of a column with linear elastic and linear or Hencky ( = natural = logarithmic ) strains
The analytical solution for these tests is:

For linear strain: ∆H = H0 * (q / E)
For logarithmic strain: ∆H = H0 * (exp (q / E) -1)

∆H : top total displacement (m)
H_0 : Original height of column (m)
q : Compressive distributed load (N/m2)
E : Young's Modulus (N/m2)
    """

    def test_same_order_elements_column_small_deformation_linear_strain(self):
        """
        Test to check if correct deformation is found using UPwSmallStrainElement2D8N and linear strain
        """

        test_name = os.path.join("test_strain_measures", "same_order_elements", "column_small_deformation_linear_strain")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve displacement top of column
        displacements = test_helper.get_displacement(simulation)

        # compare top displacement H0 q / E = 10 * -10000 / 20000 = -5 for node 3 direction 2 ( off by 1 )
        top_vertical_displacement = displacements[2][1]
        self.assertAlmostEqual( top_vertical_displacement, -5.0 )

    def test_same_order_elements_column_small_deformation_logarithmic_strain(self):
        """
        Test to check if correct deformation is found using UPwSmallStrainElement2D8N and Hencky (= natural = logarithmic) strain
        """

        test_name = os.path.join("test_strain_measures", "same_order_elements", "column_small_deformation_logarithmic_strain")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve displacement top of column
        displacements = test_helper.get_displacement(simulation)

        # compare top displacement H0 (1 - exp(q / E)) = 10 *(1 - exp(-0.5)) = -3.935 for node 3 direction 2 ( off by 1 )
        top_vertical_displacement = displacements[2][1]
        self.assertAlmostEqual( top_vertical_displacement, -3.934693402874 )

    def test_same_order_elements_column_large_deformation_logarithmic_strain(self):
        """
        Test to check if correct deformation is found using UPwUpdatedLagrangianElement2D8N and Hencky (= natural = logarithmic) strain
        """

        test_name = os.path.join("test_strain_measures", "same_order_elements", "column_large_deformation_logarithmic_strain")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve displacement top of column
        displacements = test_helper.get_displacement(simulation)

        # compare top displacement H0 (-1 + exp(q / E)) = 10 *(-1 + exp(-0.5)) = -3.935 for node 3 direction 2 ( off by 1 )
        top_vertical_displacement = displacements[2][1]
        self.assertAlmostEqual( top_vertical_displacement, -3.934693402874 )

    def test_diff_order_elements_column_small_deformation_linear_strain(self):
        """
        Test to check if correct deformation is found using SmallStrainUPwDiffOrderElement2D8N and linear strain
        """

        test_name = os.path.join("test_strain_measures", "diff_order_elements", "column_small_deformation_linear_strain")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve displacement top of column
        displacements = test_helper.get_displacement(simulation)

        # compare top displacement H0 q / E = 10 * -10000 / 20000 = -5 for node 3 direction 2 ( off by 1 )
        top_vertical_displacement = displacements[2][1]
        self.assertAlmostEqual( top_vertical_displacement, -5.0 )

    def test_diff_order_elements_column_small_deformation_logarithmic_strain(self):
        """
        Test to check if correct deformation is found using SmallStrainUPwDiffOrderElement2D8N and Hencky (= natural = logarithmic) strain
        """

        test_name = os.path.join("test_strain_measures", "diff_order_elements", "column_small_deformation_logarithmic_strain")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve displacement top of column
        displacements = test_helper.get_displacement(simulation)

        # compare top displacement H0 (1 - exp(q / E)) = 10 *(1 - exp(-0.5)) = -3.935 for node 3 direction 2 ( off by 1 )
        top_vertical_displacement = displacements[2][1]
        self.assertAlmostEqual( top_vertical_displacement, -3.934693402874 )

    def test_diff_order_elements_column_large_deformation_logarithmic_strain(self):
        """
        Test to check if correct deformation is found using UpdatedLagrangianUPwDiffOrderElement2D8N and Hencky (= natural = logarithmic) strain
        """

        test_name = os.path.join("test_strain_measures", "diff_order_elements", "column_large_deformation_logarithmic_strain")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve displacement top of column
        displacements = test_helper.get_displacement(simulation)

        # compare top displacement H0 (-1 + exp(q / E)) = 10 *(-1 + exp(-0.5)) = -3.935 for node 3 direction 2 ( off by 1 )
        top_vertical_displacement = displacements[2][1]
        self.assertAlmostEqual( top_vertical_displacement, -3.934693402874 )

if __name__ == '__main__':
    KratosUnittest.main()
