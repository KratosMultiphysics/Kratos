import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages

import test_helper

class KratosGeoMechanicsCPhiReductionProcess(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the original solution
    """
    def test_c_phi_reduction_process(self):

        # get the parameter file names for all stages
        file_path = test_helper.get_file_path('C-Phi_reduction_process')
        run_multiple_stages.run_stages(file_path, 2)

        # read results
        reader = GiDOutputFileReader()
        actual_data = reader.read_output_from(os.path.join(file_path, "stage2.post.res"))
        node_ids = [57, 151, 235]
        displacement_vectors = reader.nodal_values_at_time("DISPLACEMENT", 0.2, actual_data, node_ids)
        self.assertAlmostEqual(-0.002667, displacement_vectors[0][0])
        self.assertAlmostEqual(-0.0096777, displacement_vectors[1][0])
        self.assertAlmostEqual(-0.0115495, displacement_vectors[2][0])


if __name__ == '__main__':
    KratosUnittest.main()
