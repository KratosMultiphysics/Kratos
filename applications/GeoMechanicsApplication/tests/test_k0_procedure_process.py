import os

import KratosMultiphysics                as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsK0ProcedureProcessTests(KratosUnittest.TestCase):
    """
    This class contains tests which check the result of a K0 Procedure Process
    """

    def test_k0_procedure_k0_nc(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.6
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_nc_layers(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC,
        even when materials are stacked in layers
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc_layers")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare top layer cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.6
        sig_integrationpoint1_element6 = cauchy_stresses[5][0]
        sig_yy = sig_integrationpoint1_element6[1,1]
        sig_xx = sig_integrationpoint1_element6[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element6[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare bottom layer cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.8
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_nc_skew_layers(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC,
        even when materials are stacked in skewed layers
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc_skew_layers")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare top layer left cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.6
        sig_integrationpoint1_element3 = cauchy_stresses[2][0]
        sig_yy = sig_integrationpoint1_element3[1,1]
        sig_xx = sig_integrationpoint1_element3[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element3[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare middle layer right cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.7
        sig_integrationpoint1_element80 = cauchy_stresses[79][0]
        sig_yy = sig_integrationpoint1_element80[1,1]
        sig_xx = sig_integrationpoint1_element80[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element80[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare bottom layer left cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.8
        sig_integrationpoint1_element117 = cauchy_stresses[116][0]
        sig_yy = sig_integrationpoint1_element117[1,1]
        sig_xx = sig_integrationpoint1_element117[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element117[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_nc_skew_layers_dam(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC,
        even when materials are stacked in skewed layers
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc_skew_layers_dam")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare top layer right cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.6
        sig_integrationpoint1_element100 = cauchy_stresses[99][0]
        sig_yy = sig_integrationpoint1_element100[1,1]
        sig_xx = sig_integrationpoint1_element100[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element100[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare middle layer left cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.7
        sig_integrationpoint1_element5 = cauchy_stresses[4][0]
        sig_yy = sig_integrationpoint1_element5[1,1]
        sig_xx = sig_integrationpoint1_element5[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element5[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare bottom layer right cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.8
        sig_integrationpoint1_element51 = cauchy_stresses[50][0]
        sig_yy = sig_integrationpoint1_element51[1,1]
        sig_xx = sig_integrationpoint1_element51[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element51[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare dam1 right cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.5
        sig_integrationpoint1_element121 = cauchy_stresses[120][0]
        sig_yy = sig_integrationpoint1_element121[1,1]
        sig_xx = sig_integrationpoint1_element121[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element121[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare dam2 cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.5
        sig_integrationpoint1_element126 = cauchy_stresses[124][0]
        sig_yy = sig_integrationpoint1_element126[1,1]
        sig_xx = sig_integrationpoint1_element126[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element126[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_nc_ocr(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC and OCR
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc_ocr")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare cauchy_stress_xx = k0 * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc      = 0.6
        poisson_ur = 0.
        ocr        = 1.5
        k0 = k0_nc * ocr + ( poisson_ur / ( 1.0 - poisson_ur ) ) * ( ocr - 1.0 )
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_umat(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC = 1 - sin( PHI ),
        with PHI from UMAT material parameters
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_umat")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0. k0_nc = 1 - sin( 30 degrees )
        k0_nc = 0.5
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_simple_dike(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC = 1 - sin( PHI ) in first stage,
        with PHI from UMAT material parameters
        In second stage CAUCHY_STRESS_XX far away from the dike should remain unaltered
        Element 1562 is the right lower corner, far from the dam.
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_simple_dike")
        project_path = test_helper.get_file_path(os.path.join('.', test_name))
        cwd = os.getcwd()

        # run simulation
        n_stages = 2
        stages = test_helper.get_stages(project_path, n_stages)

        os.chdir(project_path)
        cauchy_stresses = [None] * n_stages
        for idx, stage in enumerate(stages):
            stage.Run()
            # retrieve Cauchy stress tensor of this stage
            cauchy_stresses[idx] = test_helper.get_cauchy_stress_tensor(stage)
        os.chdir(cwd)

        # compare first stage cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        # k0_nc = 1 - sin( 30 degrees )
        k0_nc = 0.5
        sig_stage1_element1562_integrationpoint1 = cauchy_stresses[0][1562-1][0]
        sig_xx_1 = sig_stage1_element1562_integrationpoint1[0,0]
        sig_yy_1 = sig_stage1_element1562_integrationpoint1[1,1]
        self.assertAlmostEqual( sig_xx_1, k0_nc*sig_yy_1 )
        sig_xy_1 = sig_stage1_element1562_integrationpoint1[0,1]
        self.assertEqual( sig_xy_1, 0.0 )

        # compare if Cauchy stress is almost unaltered far from the dam in second stage
        sig_stage2_element1562_integrationpoint1 = cauchy_stresses[1][1562-1][0]
        sig_xx_2 = sig_stage2_element1562_integrationpoint1[0,0]
        self.assertIsClose( sig_xx_2, sig_xx_1, rel_tol=0.02 )
        sig_yy_2 = sig_stage2_element1562_integrationpoint1[1,1]
        self.assertIsClose( sig_yy_2, sig_yy_1, rel_tol=0.02 )
        sig_xy_2 = sig_stage2_element1562_integrationpoint1[0,1]
        self.assertIsClose( sig_xy_2, sig_xy_1, abs_tol=1.0E01 )

if __name__ == '__main__':
    KratosUnittest.main()
