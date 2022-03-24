import sys
import os

sys.path.append(os.path.join(r'D:/kratos_new/bin/Debug'))
sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class FlowCalculations:
    def __init__(self, intrinsic_permeability=1.157e-12):
        self.intrinsic_permeability = intrinsic_permeability
        self.density_water = 1000
        self.gravity = 10
        self.viscosity = 0.001
        self.hydraulic_conductivity = self.intrinsic_permeability * self.density_water * self.gravity / self.viscosity

    def specific_dicharge(self, dphidz):
        return -1 * self.hydraulic_conductivity * dphidz

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
                output_latex_file.write(
                    f"{result_pair['value_name']} & {round(result_pair['test_result'], 2)} & "
                    f"{round(result_pair['kratos_results'], 2)} & {round(error, 2)}  \\\\ \hline \n")
        for result_pair in result_list:
            assert abs(result_pair['test_result'] - result_pair['kratos_results']) < 1e-5


class TestElementaryGroundWaterFlow(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.flow_calculations = FlowCalculations()
        self.latex_writer = LatexWriterFile()

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_hydrostatic_conditions(self):
        test_name = 'test_elementary_groundwater_flow/benchmark_1_hydrostatic.gid'
        self.latex_writer.filename = \
            test_helper.get_file_path("test_elementary_groundwater_flow/benchmark_1_hydrostatic.gid/benchmark_1_hydrostatic.tex")
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        p3 = {"value_name": "water pressure at n3 [kPa]", "test_result": -20,
              "kratos_results": test_helper.get_water_pressure(simulation)[0]}
        p4 = {"value_name": "water pressure at n4 [kPa]", "test_result": -20,
              "kratos_results": test_helper.get_water_pressure(simulation)[1]}
        phi4_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][11]
        phi3_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][8]
        phi2_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][0]
        phi1_value = test_helper.get_hydraylic_head_with_intergration_points(simulation)[2][1]
        phi4 = {"value_name": "head at n3 [m]", "test_result": 0, "kratos_results": phi4_value}
        phi3 = {"value_name": "head at n3 [m]", "test_result": 0, "kratos_results": phi3_value}
        phi2 = {"value_name": "head at n3 [m]", "test_result": 0, "kratos_results": phi2_value}
        phi1 = {"value_name": "head at n3 [m]", "test_result": 0, "kratos_results": phi1_value}
        hydraylic_disc = {"value_name": "specific discharge [m/s]", "test_result": 0,
                          "kratos_results": sum(test_helper.get_hydraulic_discharge(simulation))}
        self.latex_writer.write_latex_file_and_assert([p3, p4, phi1, phi2, phi3, phi4, hydraylic_disc])
