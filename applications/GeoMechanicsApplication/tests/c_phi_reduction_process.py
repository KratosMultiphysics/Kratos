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
        file_path = test_helper.get_file_path('C-Phi_reduction_process')
        stages = test_helper.run_stages(file_path, 2)

        # read results
        reader = test_helper.GiDOutputFileReader()
        actual_data = reader.read_output_from(os.path.join(file_path, "stage2.post.res"))
        node_ids = [57, 151, 235]
        displacement_vectors = reader.nodal_values_at_time("DISPLACEMENT", 0.2, actual_data, node_ids)
        self.assertAlmostEqual(-0.002667, displacement_vectors[0][0])
        self.assertAlmostEqual(-0.0096777, displacement_vectors[1][0])
        self.assertAlmostEqual(-0.0115495, displacement_vectors[2][0])


if __name__ == '__main__':
    KratosUnittest.main()
