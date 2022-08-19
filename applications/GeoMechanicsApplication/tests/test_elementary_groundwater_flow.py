import os
from functools import reduce
import math

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class FlowCalculations:
    """
    Analytical flow calculation class used only during testing

    """
    def __init__(self, intrinsic_permeability=1.157e-12):
        self.intrinsic_permeability = intrinsic_permeability
        self.density_water = 1000
        self.gravity = 10
        self.viscosity = 0.001
        self.hydraulic_conductivity = self.intrinsic_permeability * self.density_water * self.gravity / self.viscosity

    def specific_dicharge(self, dphidz):
        return self.hydraulic_conductivity * dphidz

    def flow_rate(self, specific_dicharge, width):
        return specific_dicharge * width / 2


class LatexWriterFile:
    """
    Class that writes results for latex documentation
    """

    def __init__(self, filename=""):
        self.filename = filename
        self.value_dict_default = {"value_name": "", "test_result": 0, "kratos_results": 0}

    def write_latex_file_and_assert(self, result_list):
        with open(self.filename, "w+") as output_latex_file:
            for result_pair in result_list:
                error = abs(result_pair['kratos_results'] - result_pair['test_result']) / \
                        (abs(result_pair['test_result']) + 1e-60)
                if result_pair["round"]:
                    output_latex_file.write(
                        f"{result_pair['value_name']} & {round(result_pair['test_result'], 2)} & "
                        f"{round(result_pair['kratos_results'], 2)} & {round(error, 2)}  \\\\ \hline \n")
                else:
                    output_latex_file.write(
                        f"{result_pair['value_name']} & {result_pair['test_result']} & "
                        f"{result_pair['kratos_results']} & {round(error, 2)}  \\\\ \hline \n")
        for result_pair in result_list:
            assert math.isclose(result_pair['test_result'], result_pair['kratos_results'], abs_tol=1e-7)


class TestElementaryGroundWaterFlow(KratosUnittest.TestCase):
    """
    Class that contains elementary groundwater flow tests.

    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.flow_calculations = FlowCalculations()
        self.latex_writer = LatexWriterFile()

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def assert_test_expectations(self, result_list):
        for result_pair in result_list:
            self.assertAlmostEqual(result_pair['test_result'], result_pair['kratos_results'], 1)

    def test_hydrostatic_conditions(self):
        """ Hydrostatic conditions test """
        test_name = 'test_elementary_groundwater_flow/benchmark_1_hydrostatic.gid'
        self.latex_writer.filename = \
            test_helper.get_file_path(
                "test_elementary_groundwater_flow/benchmark_1_hydrostatic.gid/benchmark_1_hydrostatic.tex")
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        p3 = {"value_name": "water pressure at n3 [kPa]", "test_result": -20000,
              "kratos_results": test_helper.get_water_pressure(simulation)[0], "round": True}
        p4 = {"value_name": "water pressure at n4 [kPa]", "test_result": -20000,
              "kratos_results": test_helper.get_water_pressure(simulation)[1], "round": True}
        phi4_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][11]
        phi3_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][8]
        phi2_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][0]
        phi1_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][1]
        phi4 = {"value_name": "head at n4 [m]", "test_result": 0, "kratos_results": phi4_value, "round": True}
        phi3 = {"value_name": "head at n3 [m]", "test_result": 0, "kratos_results": phi3_value, "round": True}
        phi2 = {"value_name": "head at n2 [m]", "test_result": 0, "kratos_results": phi2_value, "round": True}
        phi1 = {"value_name": "head at n1 [m]", "test_result": 0, "kratos_results": phi1_value, "round": True}
        hydraylic_disc = {"value_name": "specific discharge [m/s]", "test_result": 0,
                          "kratos_results": sum(test_helper.get_hydraulic_discharge(simulation)), "round": True}
        result_list = [p3, p4, phi1, phi2, phi3, phi4, hydraylic_disc]
        self.latex_writer.write_latex_file_and_assert(result_list)

    def test_saturated_flow_pressure_bound(self):
        """ Fully saturated soil with pressure boundary """
        test_name = 'test_elementary_groundwater_flow/benchmark_1_saturated_flow_pressure_bound.gid'
        self.latex_writer.filename = \
            test_helper.get_file_path("test_elementary_groundwater_flow/benchmark_1_saturated_flow_pressure_bound.gid/"
                                      "benchmark_1_saturated_flow_pressure_bound.tex")
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        flow_rate_1 = test_helper.get_hydraulic_discharge(simulation)[0]
        flow_rate_2 = test_helper.get_hydraulic_discharge(simulation)[1]
        flow_rate_3 = test_helper.get_hydraulic_discharge(simulation)[3]
        flow_rate_4 = test_helper.get_hydraulic_discharge(simulation)[4]
        _, y, head_list = test_helper.get_hydraylic_head_with_intergration_points(simulation)
        specific_discharge = reduce(lambda a, b: a + b if b > 0 else a, test_helper.get_hydraulic_discharge(simulation),
                                    0)

        analytical_flow_rate = self.flow_calculations.flow_rate(self.flow_calculations.specific_dicharge(1), 1)
        flow_rate_1_value = {"value_name": "flow rate Q1 [$m^{3}/s$]", "test_result": analytical_flow_rate,
                             "kratos_results": flow_rate_1, "round": False}
        flow_rate_2_value = {"value_name": "flow rate Q2 [$m^{3}/s$]", "test_result": analytical_flow_rate,
                             "kratos_results": flow_rate_2, "round": False}
        flow_rate_3_value = {"value_name": "flow rate Q3 [$m^{3}/s$]", "test_result": -1 * analytical_flow_rate,
                             "kratos_results": flow_rate_3, "round": False}
        flow_rate_4_value = {"value_name": "flow rate Q4 [$m^{3}/s$]", "test_result": -1 * analytical_flow_rate,
                             "kratos_results": flow_rate_4, "round": False}
        head_g1 = {"value_name": "hydraylic head g1 [m]", "test_result": y[3], "kratos_results": head_list[3],
                   "round": True}
        head_g2 = {"value_name": "hydraylic head g2 [m]", "test_result": y[4], "kratos_results": head_list[4],
                   "round": True}
        head_g5 = {"value_name": "hydraylic head g5 [m]", "test_result": y[0], "kratos_results": head_list[0],
                   "round": True}
        head_g6 = {"value_name": "hydraylic head g6 [m]", "test_result": y[1], "kratos_results": head_list[1],
                   "round": True}
        specific_discharge = {"value_name": "specific discharge [m/s]",
                              "test_result": self.flow_calculations.specific_dicharge(1),
                              "kratos_results": specific_discharge, "round": False}
        result_list = [flow_rate_1_value,
                       flow_rate_2_value,
                       flow_rate_3_value,
                       flow_rate_4_value,
                       head_g1,
                       head_g2,
                       head_g5,
                       head_g6,
                       specific_discharge]
        self.latex_writer.write_latex_file_and_assert(result_list)

    def test_saturated_flow_head_bound(self):
        """ Fully saturated soil with head boundary """
        test_name = 'test_elementary_groundwater_flow/benchmark_3_saturated_flow_head_bound.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        self.latex_writer.filename = \
            test_helper.get_file_path("test_elementary_groundwater_flow/benchmark_3_saturated_flow_head_bound.gid/"
                                      "benchmark_3_saturated_flow_head_bound.tex")
        simulation = test_helper.run_kratos(file_path)
        flow_rate_1 = test_helper.get_hydraulic_discharge(simulation)[0]
        flow_rate_2 = test_helper.get_hydraulic_discharge(simulation)[1]
        flow_rate_3 = test_helper.get_hydraulic_discharge(simulation)[3]
        flow_rate_4 = test_helper.get_hydraulic_discharge(simulation)[4]
        pore_pressure_1 = test_helper.get_water_pressure(simulation)[0]
        pore_pressure_2 = test_helper.get_water_pressure(simulation)[1]
        pore_pressure_3 = test_helper.get_water_pressure(simulation)[3]
        pore_pressure_4 = test_helper.get_water_pressure(simulation)[4]
        specific_discharge = reduce(lambda a, b: a + b if b > 0 else a, test_helper.get_hydraulic_discharge(simulation),
                                    0)
        analytical_flow_rate = self.flow_calculations.flow_rate(self.flow_calculations.specific_dicharge(1), 1)
        flow_rate_1_value = {"value_name": "flow rate Q1 [$m^{3}/s$]", "test_result": analytical_flow_rate,
                             "kratos_results": flow_rate_1, "round": False}
        flow_rate_2_value = {"value_name": "flow rate Q2 [$m^{3}/s$]", "test_result": analytical_flow_rate,
                             "kratos_results": flow_rate_2, "round": False}
        flow_rate_3_value = {"value_name": "flow rate Q3 [$m^{3}/s$]", "test_result": -1 * analytical_flow_rate,
                             "kratos_results": flow_rate_3, "round": False}
        flow_rate_4_value = {"value_name": "flow rate Q4 [$m^{3}/s$]", "test_result": -1 * analytical_flow_rate,
                             "kratos_results": flow_rate_4, "round": False}
        pore_pressure_1_value = {"value_name": "pore pressure p1 [Pa]", "test_result": 0,
                                 "kratos_results": pore_pressure_1, "round": True}
        pore_pressure_2_value = {"value_name": "pore pressure p2 [Pa]", "test_result": 0,
                                 "kratos_results": pore_pressure_2, "round": True}
        pore_pressure_3_value = {"value_name": "pore pressure p3 [Pa]", "test_result": 0,
                                 "kratos_results": pore_pressure_3, "round": True}
        pore_pressure_4_value = {"value_name": "pore pressure p4 [Pa]", "test_result": 0,
                                 "kratos_results": pore_pressure_4, "round": True}
        specific_discharge = {"value_name": "specific discharge [m/s]",
                              "test_result": self.flow_calculations.specific_dicharge(1),
                              "kratos_results": specific_discharge, "round": False}
        result_list = [flow_rate_1_value,
                       flow_rate_2_value,
                       flow_rate_3_value,
                       flow_rate_4_value,
                       pore_pressure_1_value,
                       pore_pressure_2_value,
                       pore_pressure_3_value,
                       pore_pressure_4_value,
                       specific_discharge]
        self.latex_writer.write_latex_file_and_assert(result_list)

    def test_saturated_flux_bound(self):
        """ Fully saturated soil with flux boundary """
        test_name = 'test_elementary_groundwater_flow/benchmark_4_saturated_flux_bound.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        flow_rate_1 = test_helper.get_hydraulic_discharge(simulation)[0]
        flow_rate_2 = test_helper.get_hydraulic_discharge(simulation)[1]
        flow_rate_3 = test_helper.get_hydraulic_discharge(simulation)[3]
        flow_rate_4 = test_helper.get_hydraulic_discharge(simulation)[4]
        self.assertTrue(math.isclose(abs(flow_rate_1), abs(self.flow_calculations.flow_rate(self.flow_calculations.specific_dicharge(1), 1))))
        self.assertTrue(math.isclose(abs(flow_rate_2), abs(self.flow_calculations.flow_rate(self.flow_calculations.specific_dicharge(1), 1))))
        self.assertTrue(math.isclose(abs(flow_rate_3), abs(self.flow_calculations.flow_rate(self.flow_calculations.specific_dicharge(1), 1))))
        self.assertTrue(math.isclose(abs(flow_rate_4), abs(self.flow_calculations.flow_rate(self.flow_calculations.specific_dicharge(1), 1))))
        self.assertTrue(math.isclose(test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][8], -1.5))
        self.assertTrue(math.isclose(test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][11], -1.5))