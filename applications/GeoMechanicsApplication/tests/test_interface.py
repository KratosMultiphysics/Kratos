import sys
import os

import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsInterfaceTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which check if interfaces between soil and structural elements are
    calculated correctly
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_interface_side_cohesive(self):
        """
        Tests an interface between a fixed line-element with and a non-fixed surface.

        The surface has a prescribed vertical displacement of -0.1 m and a horizontal line load of -1 kN
        The interface has a cohesion of 10 kN.
        Expected shear stress in the interface is 10 kN

        :return:
        """
        test_name = 'test_interface_side_cohesive'
        file_path = test_helper.get_file_path(os.path.join('.','test_interfaces', test_name + '.gid'))

        simulation = test_helper.run_kratos(file_path)
        local_stress_vector = test_helper.get_local_stress_vector(simulation)
        local_stress_vector_interface = [element for element in local_stress_vector if len(element) == 4]

        x_stress = 10000
        precision = 0.01 * x_stress  # 1% precision

        for local_stress_quad in local_stress_vector_interface:
            error_gauss_1 = abs(x_stress - local_stress_quad[0][0])
            error_gauss_2 = abs(x_stress - local_stress_quad[1][0])

            self.assertLess(error_gauss_1, precision)
            self.assertLess(error_gauss_2, precision)

    def test_interface_on_beam(self):
        """
        Tests an interface on a beam. In this test a calculation is done with and without interface,
        where the interface has a cohesion of 1000 kN and a stiffness of 1e12 Pa. Results should be approximately equal

        :return:
        """

        # calculate reference case without interface
        test_name = 'test_beam'
        file_path = test_helper.get_file_path(os.path.join('.','test_interfaces', test_name + '.gid'))

        simulation = test_helper.run_kratos(file_path)
        base_displacement = test_helper.get_displacement(simulation)

        base_x_displacement = [displacement[0] for displacement in base_displacement]
        base_y_displacement = [displacement[1] for displacement in base_displacement]

        max_base_x_displacement = max(base_x_displacement)
        max_base_y_displacement, min_base_y_displacement = max(base_y_displacement), min(base_y_displacement)

        # calculate case with strong interface
        test_name = 'test_interface_on_beam'
        file_path = test_helper.get_file_path(
            os.path.join('.', 'test_interfaces', test_name + '.gid'))

        simulation = test_helper.run_kratos(file_path)
        interface_displacement = test_helper.get_displacement(simulation)

        interface_x_displacement = [displacement[0] for displacement in interface_displacement]
        interface_y_displacement = [displacement[1] for displacement in interface_displacement]

        max_interface_x_displacement = max(interface_x_displacement)
        max_interface_y_displacement, min_interface_y_displacement = max(interface_y_displacement), min(interface_y_displacement)

        # assert base displacement is approximately equal to displacement with interface
        precision = 0.01
        self.assertLess(abs(max_interface_x_displacement - max_base_x_displacement), precision*max_base_x_displacement)
        self.assertLess(abs(max_interface_y_displacement - max_base_y_displacement), precision*max_base_y_displacement)
        self.assertLess(abs(min_interface_y_displacement - min_base_y_displacement), precision*max_base_y_displacement)

    def test_weak_interface_on_beam(self):
        """
        Tests an interface on a beam. In this test a calculation is done with a very weak interface.
        Displacement in the beam element should be much greater than the displacement in the soil.

        :return:
        """

        # calculate case with strong interface
        test_name = 'test_weak_interface_on_beam'
        file_path = test_helper.get_file_path(
            os.path.join('.', 'test_interfaces', test_name + '.gid'))

        simulation = test_helper.run_kratos(file_path)
        model_part = simulation._list_of_output_processes[0].model_part

        # Get beam and soil elements
        beam_elements = [element for element in model_part.Elements if element.GetGeometry().PointsNumber() == 2]
        soil_elements = [element for element in model_part.Elements if element.GetGeometry().PointsNumber() == 3]

        # Get average x-displacement beam- and soil- elements
        x_displacement_beam_elements = [test_helper.compute_mean_list([node.GetSolutionStepValue(Kratos.DISPLACEMENT)[0]
                                         for node in beam_element.GetNodes()]) for beam_element in beam_elements]
        x_displacement_soil_elements = [test_helper.compute_mean_list([node.GetSolutionStepValue(Kratos.DISPLACEMENT)[0]
                                         for node in soil_element.GetNodes()]) for soil_element in soil_elements]

        # Get max x-displacement beam- and soil- elements
        max_x_disp_beam = max(x_displacement_beam_elements)
        max_x_disp_soil = max(x_displacement_soil_elements)

        # Assert if beam displacement is >> than soil displacement
        self.assertGreater(max_x_disp_beam, 1e8 * max_x_disp_soil)

if __name__ == '__main__':
    KratosUnittest.main()
