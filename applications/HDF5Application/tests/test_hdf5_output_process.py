from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5

import hdf5_output_process

class TestHDF5OutputProcess(KratosUnittest.TestCase):
    def setUp(self):
        pass

     def test_HDF5SerialOutput(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
