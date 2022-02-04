import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))
sys.path.append(os.path.join('D:\\kratos'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class TestPipingElements(KratosUnittest.TestCase):
    """
    This class contains test that verifies if the piping elements work
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_different_heights_pipe_elements(self):
        tests_name = ['test_4_1', 'test_4_2', 'test_4_3', 'test_4_4', 'test_4_5']
        for test_name in tests_name:
            file_path = test_helper.get_file_path(
                os.path.join('./SteadyStatePipeElementWithEmbankment/test_different_heights_pipe_elements',
                             test_name + '.gid'))
            simulation = test_helper.run_kratos(file_path)
            x, y, head = test_helper.get_hydraylic_head_with_intergration_points(simulation)
            csv_path = os.path.join(file_path, test_name + '.csv')
            with open(csv_path, 'w', encoding='UTF8', newline='') as f:
                f.write('x,y,head\n')
                for counter, value in enumerate(head):
                    f.write(str(x[counter]) + ',' + str(y[counter]) + ',' + str(value) + '\n')

    def test_SteadyStatePipeElementWithEmbankment(self):
        tests_name = ['SteadyStatePipeElementWithEmbankment_4_1',
                      'SteadyStatePipeElementWithEmbankment_4_2',
                      'SteadyStatePipeElementWithEmbankment_4_3',
                      'SteadyStatePipeElementWithEmbankment_4_4',
                      'SteadyStatePipeElementWithEmbankment_4_5',
                      'SteadyStatePipeElementWithEmbankment_4_6']

        for test_name in tests_name:
            file_path = test_helper.get_file_path(
                os.path.join('./SteadyStatePipeElementWithEmbankment', test_name + '.gid'))
            simulation = test_helper.run_kratos(file_path)
            model_part = simulation._list_of_output_processes[0].model_part


if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small']  # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPipingElements]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
