import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

class TestNDDataIO(UnitTest.TestCase):
    def setUp(self) -> None:
        self.data_communicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
        if self.data_communicator.IsDistributed():
            parameters = Kratos.Parameters("""{
                "file_name" : "test_nd_data_io.h5",
                "file_access_mode": "truncate",
                "file_driver": "mpio",
                "echo_level" : 0
            }""")
        else:
            parameters = Kratos.Parameters("""{
                "file_name" : "test_nd_data_io.h5",
                "file_access_mode": "truncate",
                "file_driver": "sec2",
                "echo_level" : 0
            }""")
        self.h5_file = KratosHDF5.HDF5File(self.data_communicator, parameters)

    def tearDown(self) -> None:
        self.data_communicator.Barrier()
        DeleteFileIfExisting("test_nd_data_io.h5")

    def test_BoolNDData(self):
        numpy_array = np.random.random((self.data_communicator.Rank() + 10, 5, 6)) > 0.5
        write_nd_data = Kratos.BoolNDData(numpy_array, copy=False)
        write_attributes = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")

        nd_data_io = KratosHDF5.NDDataIO(Kratos.Parameters("""{"prefix": "/nd_data/"}"""), self.h5_file)
        nd_data_io.Write("bool_test", write_nd_data, write_attributes)

        read_nd_data, read_attributes = nd_data_io.Read("bool_test")

        self.assertTrue(isinstance(read_nd_data, Kratos.UIntNDData)) # bools are always read as uint from HDF5
        self.assertTrue(read_attributes.IsEquivalentTo(write_attributes))
        self.assertTrue(np.all(read_nd_data.to_numpy() == numpy_array))

    def test_UIntNDData(self):
        numpy_array = (np.random.random((self.data_communicator.Rank() + 10, 5, 6)) * 255.0).astype(np.uint8)
        write_nd_data = Kratos.UIntNDData(numpy_array, copy=False)
        write_attributes = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")

        nd_data_io = KratosHDF5.NDDataIO(Kratos.Parameters("""{"prefix": "/nd_data/"}"""), self.h5_file)
        nd_data_io.Write("uint_test", write_nd_data, write_attributes)

        read_nd_data, read_attributes = nd_data_io.Read("uint_test")

        self.assertTrue(isinstance(read_nd_data, Kratos.UIntNDData))
        self.assertTrue(read_attributes.IsEquivalentTo(write_attributes))
        self.assertTrue(np.all(read_nd_data.to_numpy() == numpy_array))

    def test_IntNDData(self):
        numpy_array = (np.random.random((self.data_communicator.Rank() + 10, 5, 6)) * 1000).astype(np.int32)
        write_nd_data = Kratos.IntNDData(numpy_array, copy=False)
        write_attributes = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")

        nd_data_io = KratosHDF5.NDDataIO(Kratos.Parameters("""{"prefix": "/nd_data/"}"""), self.h5_file)
        nd_data_io.Write("uint_test", write_nd_data, write_attributes)

        read_nd_data, read_attributes = nd_data_io.Read("uint_test")

        self.assertTrue(isinstance(read_nd_data, Kratos.IntNDData))
        self.assertTrue(read_attributes.IsEquivalentTo(write_attributes))
        self.assertTrue(np.all(read_nd_data.to_numpy() == numpy_array))

    def test_FloatNDData(self):
        numpy_array = (np.random.random((self.data_communicator.Rank() + 10, 5, 6)) * 1000).astype(np.float64)
        write_nd_data = Kratos.DoubleNDData(numpy_array, copy=False)
        write_attributes = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")

        nd_data_io = KratosHDF5.NDDataIO(Kratos.Parameters("""{"prefix": "/nd_data/"}"""), self.h5_file)
        nd_data_io.Write("uint_test", write_nd_data, write_attributes)

        read_nd_data, read_attributes = nd_data_io.Read("uint_test")

        self.assertTrue(isinstance(read_nd_data, Kratos.DoubleNDData))
        self.assertTrue(read_attributes.IsEquivalentTo(write_attributes))
        self.assertTrue(np.all(read_nd_data.to_numpy() == numpy_array))

    def test_EmptyRankNDData(self):
        numpy_array = (np.random.random(((self.data_communicator.Rank() + 1) % 3, 5, 6)) * 1000).astype(np.float64)
        write_nd_data = Kratos.DoubleNDData(numpy_array, copy=False)
        write_attributes = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")

        nd_data_io = KratosHDF5.NDDataIO(Kratos.Parameters("""{"prefix": "/nd_data/"}"""), self.h5_file)
        nd_data_io.Write("uint_test", write_nd_data, write_attributes)

        read_nd_data, read_attributes = nd_data_io.Read("uint_test")

        self.assertTrue(isinstance(read_nd_data, Kratos.DoubleNDData))
        self.assertTrue(read_attributes.IsEquivalentTo(write_attributes))
        self.assertTrue(np.all(read_nd_data.to_numpy() == numpy_array))

if __name__ == "__main__":
    UnitTest.main()