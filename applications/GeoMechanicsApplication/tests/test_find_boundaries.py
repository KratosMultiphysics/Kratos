import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsMasterSlaveConstraints(KratosUnittest.TestCase):

    def test_find_boundaries(self):
        test_files_path = test_helper.get_file_path("find_boundaries")
        test_helper.run_kratos(test_files_path)

if __name__ == "__main__":
    KratosUnittest.main()
