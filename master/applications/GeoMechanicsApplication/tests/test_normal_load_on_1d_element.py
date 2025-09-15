import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosStruct
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsNormalLoad1DTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests to test if normal loads are correctly calculated on 1D elements.
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def calculate_analytical_deflection_simply_supported_beam_with_normal_load(self, properties, L, n_elem, load):
        """
        Calculates the deflection of a simply supported beam following a normal load
        Parameters
        ----------
        properties: Kratos Properties
        L:  Length Beam
        n_elem: number of elements beam
        load: normal load value

        Returns
        -------

        """

        # get parameters
        E = properties[Kratos.YOUNG_MODULUS]
        I = properties[KratosStruct.I22]

        # discretise x coordinate
        coords = [i * L/n_elem for i in range(n_elem + 1)]

        # calculate displacement at each coordinate
        displacements = [load * x / (24 * E * I) * (L ** 3 - 2 * L * x ** 2 + x ** 3) for x in coords]
        return displacements

    def test_normal_load_on_beam(self):
        """
        Tests the deflection of a beam following a normal load.
        Returns
        -------

        """
        # read test
        test_name = 'test_normal_load_on_beam'
        file_path = test_helper.get_file_path(os.path.join(test_name + '.gid'))

        # run calculation
        simulation = test_helper.run_kratos(file_path)

        # get results
        nodes = simulation.model.GetModelPart("PorousDomain.beam").Nodes
        y_displacements = [node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) for node in nodes]

        # calculate expected result
        properties = simulation.model.GetModelPart("PorousDomain.beam").Properties[1]
        expected_result = self.calculate_analytical_deflection_simply_supported_beam_with_normal_load\
            (properties, 3, len(y_displacements)-1, -10)

        # calculate relative difference
        relative_difference = [(y_displacement - expected) / expected for y_displacement, expected
                               in zip(y_displacements, expected_result) if expected != 0]

        # assert relative differences
        for diff in relative_difference:
            self.assertAlmostEqual(diff, 0, places=3)

    @KratosUnittest.skip("unit test skipped due to problem in the core, remove skip statement after #10950 is solved")
    def test_normal_load_on_beam_higher_order(self):
        """
        Tests the deflection of a beam following a normal load.
        Returns
        -------

        """
        # read test
        test_name = 'test_normal_load_on_beam_higher_order'
        file_path = test_helper.get_file_path(os.path.join(test_name + '.gid'))

        # run calculation
        simulation = test_helper.run_kratos(file_path)

        # get results
        nodes = simulation.model.GetModelPart("PorousDomain.beam").Nodes
        y_displacements = [node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) for node in nodes]

        # calculate expected result
        properties = simulation.model.GetModelPart("PorousDomain.beam").Properties[1]
        expected_result = self.calculate_analytical_deflection_simply_supported_beam_with_normal_load\
            (properties, 3, len(y_displacements)-1, -10)

        # calculate relative difference
        relative_difference = [(y_displacement - expected) / expected for y_displacement, expected
                               in zip(y_displacements, expected_result) if expected != 0]

        # assert relative differences
        for diff in relative_difference:
            self.assertAlmostEqual(diff, 0, places=3)

    def test_normal_load_on_beam_and_soil(self):
        """
        Tests the deflection of a beam following a normal load. Where the beam is attached to soil with a very low
        stiffness
        Returns
        -------

        """

        # read test
        test_name = 'test_normal_load_on_beam_and_soil'
        file_path = test_helper.get_file_path(os.path.join(test_name + '.gid'))

        # run calculation
        simulation = test_helper.run_kratos(file_path)

        # get results
        nodes = simulation.model.GetModelPart("PorousDomain.beam").Nodes
        y_displacements = [node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) for node in nodes]

        # calculate expected result
        properties = simulation.model.GetModelPart("PorousDomain.beam").Properties[2]
        expected_result = self.calculate_analytical_deflection_simply_supported_beam_with_normal_load\
            (properties, 3, len(y_displacements)-1, -10)

        # calculate relative difference
        relative_difference = [(y_displacement - expected) / expected for y_displacement, expected
                               in zip(y_displacements, expected_result) if expected != 0]

        # assert relative differences
        for diff in relative_difference:
            self.assertAlmostEqual(diff, 0, places=3)


if __name__ == '__main__':
    KratosUnittest.main()
