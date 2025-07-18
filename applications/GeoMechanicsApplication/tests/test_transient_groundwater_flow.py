import sys
import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsTransientGroundWaterFlowTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """
    def test_Transient_Case_B1_2D3N(self):
        test_name = 'test_Transient_Case_B1_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        stages = run_multiple_stages.run_stages(file_path, 2)

        water_pressure_stage_2 = test_helper.get_water_pressure(stages[-1])

        p_height_1_72_m = water_pressure_stage_2[32]

        self.assertAlmostEqual(1.8902712667818171, p_height_1_72_m, 3)

    def test_Transient_Case_A1_2D3N(self):
        test_name = 'test_Transient_Case_A1_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        stages = run_multiple_stages.run_stages(file_path, 2)

        water_pressure_stage_2 = test_helper.get_water_pressure(stages[-1])

        p_height_1_72_m = water_pressure_stage_2[32]

        self.assertAlmostEqual(3.730389835521539, p_height_1_72_m)

    def test_Transient_Case_A1_2D6N(self):
        test_name = 'test_Transient_Case_A1_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        stages = run_multiple_stages.run_stages(file_path, 2)

        water_pressure = test_helper.get_water_pressure(stages[-1])

        p_height_1_72_m = water_pressure[87]

        self.assertAlmostEqual(2.617220832597392, p_height_1_72_m)


if __name__ == '__main__':
    KratosUnittest.main()
