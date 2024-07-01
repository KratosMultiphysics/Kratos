# from KratosMultiphysics import * as Kratos

import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

import test_helper

class KratosGeoMechanicsCPhiReductionProcess(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the original solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass
        
    def test_c_phi_reduction_process(self):

        # get the parameter file names for all stages
        test_name = 'C-Phi_reduction_process'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        stages = test_helper.run_stages(file_path, 2)

        # read results
        reader = test_helper.GiDOutputFileReader()
        reader.read_output_from(os.path.join(file_path, "stage2.post.res"))
        displacement = test_helper.get_displacement(stages[1])
        self.assertAlmostEqual(-0.0008964807371535437, displacement[12][0])
        self.assertAlmostEqual(-0.0011135191709369784, displacement[27][0])
        self.assertAlmostEqual(-0.00967769824048078, displacement[150][0])


if __name__ == '__main__':
    KratosUnittest.main()
