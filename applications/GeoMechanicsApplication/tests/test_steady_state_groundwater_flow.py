import sys
import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsSteadyStateGroundWaterFlowTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.test_confined_aquifer = [('test_flow_under_dam/test_pressure_in_confined_aquifer_higher', 5.499),
                                      ('test_flow_under_dam/test_pressure_in_confined_aquifer_lower', 0.423),
                                      ('test_flow_under_dam/test_pressure_in_confined_aquifer_reversed', 2.54028)]


    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_SteadyState_OneDFlow_2D3N(self):
        test_name = 'test_SteadyState_OneDFlow_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 1.0)

    def test_SteadyState_DamConfinedFlow_2D3N(self):
        test_name = 'test_SteadyState_DamConfinedFlow_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 1.11309991)

    def test_SteadyState_DamConfinedFlow_2D6N(self):
        test_name = 'test_SteadyState_DamConfinedFlow_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 1.08514372)

    def test_SteadyState_DamConfinedFlowWithSheetPile_2D6N(self):
        test_name = 'test_SteadyState_DamConfinedFlowWithSheetPile_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        self.assert_outflow_discharge(simulation, 0.827800245)
    
    def test_darcy_law_on_one_element(self):
        test_name = 'test_darcy_law'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        # analytical result
        Q = 2
        # results from kratos
        outflow_discharge = self.calculate_outflow_discharge(simulation)
        self.assert_outflow_discharge(simulation, Q)
        error_outflow_discharge = abs(outflow_discharge - Q) / (Q + 1e-60)
        print('Writing tex file in: ', os.path.abspath(file_path + "\\test_darcy_law_on_one_element.tex"))
        output_file_for_latex = open(file_path + "\\test_darcy_law_on_one_element.tex", "w")
        output_file_for_latex.write(' & '.join(['Q',
                                                str(round(Q, 2)),
                                                str(round(outflow_discharge, 2)),
                                                str(round(error_outflow_discharge * 100, 2))]) +
                                    ' \\\\ \hline \n')
        output_file_for_latex.close()

    def test_flow_under_dam(self):
        for test_name, Q in self.test_confined_aquifer:
            with self.subTest():
                file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
                simulation = test_helper.run_kratos(file_path)
                outflow_discharge = self.calculate_outflow_discharge(simulation)
                error_outflow_discharge = abs(outflow_discharge - Q) / (Q + 1e-60)
                self.assertTrue(error_outflow_discharge < 0.03)
                print('Writing tex file in: ', os.path.abspath(file_path + "\\test_flow_under_dam.tex"))
                output_file_for_latex = open(file_path + "\\test_flow_under_dam.tex", "w")
                output_file_for_latex.write(' & '.join(['Q',
                                                        str(round(Q, 2)),
                                                        str(round(outflow_discharge, 2)),
                                                        str(round(error_outflow_discharge * 100, 2))]) +
                                            ' \\\\ \hline \n')
                output_file_for_latex.close()

    def test_flow_rate_heterogeneous_soil(self):
        test_name = 'flow_rate_heterogeneous_soil'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        outflow_discharge = self.calculate_outflow_discharge(simulation)
        analytical_solution_outflow_discharge = 0.056
        error_outflow_discharge = abs(outflow_discharge - analytical_solution_outflow_discharge) / \
                                  (analytical_solution_outflow_discharge + 1e-60)
        self.assertTrue(error_outflow_discharge < 0.03)
        print('Writing tex file in: ', os.path.abspath(file_path + "\\test_pressure_in_confined_aquifer.tex"))
        output_file_for_latex = open(file_path + "\\test_flow_rate_heterogeneous_soil.tex","w")
        output_file_for_latex.write(' & '.join(['Q',
                                                str(round(analytical_solution_outflow_discharge,2)),
                                                str(round(outflow_discharge,2)),
                                                str(round(error_outflow_discharge*100,2))]) +
                                                ' \\\\ \hline \n')
        output_file_for_latex.close()

    def test_pressure_in_confined_aquifer(self):
        # run kratos model
        test_name = 'test_pressure_in_confined_aquifer'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        # run test for outflow discharge
        analytical_solution_outflow_discharge = 2.54028
        outflow_discharge = self.calculate_outflow_discharge(simulation)
        error_outflow_discharge = abs(outflow_discharge - analytical_solution_outflow_discharge) / \
                                  (analytical_solution_outflow_discharge + 1e-60)
        self.assertTrue(error_outflow_discharge < 0.03)
        analytical_solution_water_pressure = -127.53
        # get water pressure
        water_pressure = test_helper.get_water_pressure(simulation)
        pore_pressure_middle = water_pressure[343]
        error_water_pressure = abs(pore_pressure_middle - analytical_solution_water_pressure) / \
                                  (abs(analytical_solution_water_pressure) + 1e-60)
        self.assertTrue(error_outflow_discharge < 0.03)
        print('Writing tex file in: ', os.path.abspath(file_path + "\\test_pressure_in_confined_aquifer.tex"))
        output_file_for_latex = open(file_path + "\\test_pressure_in_confined_aquifer.tex","w")
        output_file_for_latex.write(' & '.join(['water pressure',
                                                str(round(analytical_solution_water_pressure,2)),
                                                str(round(pore_pressure_middle,2)),
                                                str(round(error_water_pressure*100,2))]) + ' \\\\ \hline \n')
        output_file_for_latex.write(' & '.join(['Q',
                                                str(round(analytical_solution_outflow_discharge,2)),
                                                str(round(outflow_discharge,2)),
                                                str(round(error_outflow_discharge*100,2))]) +
                                                ' \\\\ \hline \n')
        output_file_for_latex.close()

    def calculate_outflow_discharge(self, simulation):
        """
        calculate total outflow discharge.

        :param simulation: Kratos simulation
        :return: total outflow discharge
        """

        hydraulic_discharge = test_helper.get_hydraulic_discharge(simulation)

        from functools import reduce
        return reduce(lambda a,b: a+b if b>0 else a, hydraulic_discharge, 0)

    def calculate_total_discharge(self, simulation):
        """
        calculate total outflow discharge.

        :param simulation: Kratos simulation
        :return: total outflow discharge
        """

        hydraulic_discharge = test_helper.get_hydraulic_discharge(simulation)

        from functools import reduce
        return reduce(lambda a,b: a+b, hydraulic_discharge, 0)


    def assert_outflow_discharge(self, simulation, expected_value):
        """
        assert total outflow discharge.

        :param simulation: Kratos simulation
        :param expected_value: expected value
        :return:
        """

        outflow_discharge = self.calculate_outflow_discharge(simulation)
        error = abs(outflow_discharge - expected_value) / (expected_value + 1e-60)
        tolerated_error = 1e-6
        self.assertLess(error, tolerated_error)

    def assert_total_discharge(self, simulation, max_value):
        """
        assert total outflow discharge.

        :param simulation: Kratos simulation
        :param expected_value: expected value
        :return:
        """

        total_discharge = abs(self.calculate_total_discharge(simulation))
        self.assertLess(total_discharge, max_value)


if __name__ == '__main__':
    KratosUnittest.main()
