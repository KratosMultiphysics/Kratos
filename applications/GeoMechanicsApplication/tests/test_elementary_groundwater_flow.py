import sys
import os


import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class FlowCalculations:
    def __init__(self, intrinsic_permeability=1.157e-12):
        self.intrinsic_permeability = intrinsic_permeability
        self.density_water = 1000
        self.gravity = 10
        self.viscosity = 0.001
        self.hydraulic_conductivity =self.intrinsic_permeability * self.density_water * self.gravity / self.viscosity

    def specific_dicharge(self, dphidz):
        return -1 * self.hydraulic_conductivity * dphidz

    def flow_rate(self, specific_dicharge, width):
        return specific_dicharge * width / 2

class TestElementaryGroundWaterFlow(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.flow_calculations = FlowCalculations()

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_hydrostatic_conditions(self):
        test_name = 'test_elementary_groundwater_flow/benchmark_1_hydrostatic.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        p4 = test_helper.get_water_pressure(simulation)[0]
        p5 = test_helper.get_water_pressure(simulation)[1]
        assert abs(p4 - (-20)) < 1e-5
        assert abs(p5 - (-20)) < 1e-5
        assert abs(test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][11] - 0) < 1e-5
        assert abs(test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][8] - 0) < 1e-5
        assert abs(test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][0] - 0) < 1e-5
        assert abs(test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][1] - 0) < 1e-5
        assert sum(test_helper.get_hydraulic_discharge(simulation)) == 0
