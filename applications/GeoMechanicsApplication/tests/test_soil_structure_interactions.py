import sys
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsSoilStructureInteractionTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass


    def test_truss_between_soil(self):
        """
        Two blocks of soil are attached with very stiff truss. A point load is applied to 1 soil block.
        Check if displacements at both sides of the truss are equal.

        :return:
        """
        test_name = 'test_truss_between_soils'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)

        # check if displacement in soil at both sides of truss is equal
        self.assertAlmostEqual(displacements[6][0], displacements[8][0])

    def test_cable_between_soil(self):
        """
        Two blocks of soil are attached with very stiff cable. A point load is applied to 1 soil block.
        Check if displacements at both sides of the cable are equal.

        :return:
        """

        test_name = 'test_cable_between_soils'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)

        # check if displacement in soil at both sides of truss is equal
        self.assertAlmostEqual(displacements[6][0], displacements[8][0])

if __name__ == '__main__':
    KratosUnittest.main()
