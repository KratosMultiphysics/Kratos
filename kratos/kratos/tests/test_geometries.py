from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

class TestModelPart(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def test_tetrahedra_3D4N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D4N( ) )

    def test_tetrahedra_2D3N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D3N( ) )

    def test_tetrahedra_2D6N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTriangle2D6N( ) )

    def test_tetrahedra_3D4N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D4N( ) )

    def test_tetrahedra_3D10N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestTetrahedra3D10N( ) )

    def test_tetrahedra_3D8N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D8N( ) )

    def test_tetrahedra_3D27N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D27N( ) )
        
    def test_tetrahedra_3D20N(self):
        tester = GeometryTesterUtility()
        self.assertTrue( tester.TestHexahedra3D20N( ) )
        


if __name__ == '__main__':
    KratosUnittest.main()
