import numpy as np
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

class TestTensorAdaptorIO(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        ReadModelPart(str(TestTensorAdaptorIO.GetInputMDPAPath()), cls.model_part)

        def set_values(container):
            for entity in container:
                entity.SetValue(Kratos.PRESSURE, entity.Id)
                entity.SetValue(Kratos.VELOCITY, Kratos.Array3([entity.Id, entity.Id + 1, entity.Id + 2]))
                entity.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, Kratos.Vector(4, entity.Id))
                entity.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix(2, 3, entity.Id))


        set_values(cls.model_part.Nodes)
        set_values(cls.model_part.Conditions)
        set_values(cls.model_part.Elements)

        cls.data_communicator: Kratos.DataCommunicator = cls.model_part.GetCommunicator().GetDataCommunicator()

    def setUp(self) -> None:
        if self.data_communicator.IsDistributed():
            parameters = Kratos.Parameters("""{
                "file_name" : "test_tensor_adaptor_io.h5",
                "file_access_mode": "truncate",
                "file_driver": "mpio",
                "echo_level" : 0
            }""")
        else:
            parameters = Kratos.Parameters("""{
                "file_name" : "test_tensor_adaptor_io.h5",
                "file_access_mode": "truncate",
                "file_driver": "sec2",
                "echo_level" : 0
            }""")
        self.h5_file = KratosHDF5.HDF5File(self.data_communicator, parameters)

    def tearDown(self) -> None:
        self.data_communicator.Barrier()
        DeleteFileIfExisting("test_tensor_adaptor_io.h5")

    def __TestTensorAdaptor(self, variable, container):
        attributes_in = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(container, variable)
        ta.Check()
        ta.CollectData()

        tensor_adaptor_io = KratosHDF5.TensorAdaptorIO(Kratos.Parameters("""{"prefix": "/tensor_adaptors/"}"""), self.h5_file)
        tensor_adaptor_io.Write(f"{container.__class__.__name__}", ta, attributes_in)

        ta_out = Kratos.TensorAdaptors.VariableTensorAdaptor(container, variable)
        ta_attributes_out = tensor_adaptor_io.Read(f"{container.__class__.__name__}", ta_out)
        self.assertTrue(ta_attributes_out.IsEquivalentTo(attributes_in))
        self.assertEqual(np.linalg.norm(ta_out.data - ta.data), 0)

    def __TestWriteVariableAndReadTensorAdaptor(self, variable, hdf5_io_type, container, prefix: str):
        params = Kratos.Parameters("""{"prefix": ""}""")
        params["prefix"].SetString(prefix)
        params.AddStringArray("list_of_variables", [variable.Name()])
        io = hdf5_io_type(params, self.h5_file)
        io.Write(self.model_part)

        ta_orig = Kratos.TensorAdaptors.VariableTensorAdaptor(container, variable)
        ta_orig.Check()
        ta_orig.CollectData()

        tensor_io_params = Kratos.Parameters("""{"prefix": ""}""")
        tensor_io_params["prefix"].SetString(prefix)
        ta_read = Kratos.TensorAdaptors.VariableTensorAdaptor(container, variable)
        _ = KratosHDF5.TensorAdaptorIO(tensor_io_params, self.h5_file).Read(variable.Name(), ta_read)

        # following is required since we flatten the higher dimensions than 1st dimension in the data shape when we write the variables.
        orig_data_shape = ta_orig.DataShape()
        current_shape = [ta_read.Shape()[0]]
        current_shape.extend(orig_data_shape)

        self.assertEqual(np.linalg.norm(ta_orig.data - ta_read.data.reshape(current_shape)), 0)

    def test_ReadWriteTensorAdaptorScalar(self):
        self.__TestTensorAdaptor(Kratos.PRESSURE, self.model_part.Nodes)
        self.__TestTensorAdaptor(Kratos.PRESSURE, self.model_part.Conditions)
        self.__TestTensorAdaptor(Kratos.PRESSURE, self.model_part.Elements)

    def test_ReadWriteTensorAdaptorArray3(self):
        self.__TestTensorAdaptor(Kratos.VELOCITY, self.model_part.Nodes)
        self.__TestTensorAdaptor(Kratos.VELOCITY, self.model_part.Conditions)
        self.__TestTensorAdaptor(Kratos.VELOCITY, self.model_part.Elements)

    def test_ReadWriteTensorAdaptorVector(self):
        self.__TestTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, self.model_part.Nodes)
        self.__TestTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, self.model_part.Conditions)
        self.__TestTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, self.model_part.Elements)

    def test_ReadWriteTensorAdaptorMatrix(self):
        self.__TestTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part.Nodes)
        self.__TestTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part.Conditions)
        self.__TestTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part.Elements)

    def test_WriteVariableReadTensorAdaptorScalar(self):
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.PRESSURE, KratosHDF5.HDF5NodalDataValueIO, self.model_part.GetCommunicator().LocalMesh().Nodes, "/ResultsData/NodalDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.PRESSURE, KratosHDF5.HDF5ConditionDataValueIO, self.model_part.GetCommunicator().LocalMesh().Conditions, "/ResultsData/ConditionDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.PRESSURE, KratosHDF5.HDF5ElementDataValueIO, self.model_part.GetCommunicator().LocalMesh().Elements, "/ResultsData/ElementDataValues/")

    def test_WriteVariableReadTensorAdaptorArray3(self):
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.VELOCITY, KratosHDF5.HDF5NodalDataValueIO, self.model_part.GetCommunicator().LocalMesh().Nodes, "/ResultsData/NodalDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.VELOCITY, KratosHDF5.HDF5ConditionDataValueIO, self.model_part.GetCommunicator().LocalMesh().Conditions, "/ResultsData/ConditionDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.VELOCITY, KratosHDF5.HDF5ElementDataValueIO, self.model_part.GetCommunicator().LocalMesh().Elements, "/ResultsData/ElementDataValues/")

    def test_WriteVariableReadTensorAdaptorVector(self):
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, KratosHDF5.HDF5NodalDataValueIO, self.model_part.GetCommunicator().LocalMesh().Nodes, "/ResultsData/NodalDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, KratosHDF5.HDF5ConditionDataValueIO, self.model_part.GetCommunicator().LocalMesh().Conditions, "/ResultsData/ConditionDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, KratosHDF5.HDF5ElementDataValueIO, self.model_part.GetCommunicator().LocalMesh().Elements, "/ResultsData/ElementDataValues/")

    def test_WriteVariableReadTensorAdaptorMatrix(self):
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, KratosHDF5.HDF5NodalDataValueIO, self.model_part.GetCommunicator().LocalMesh().Nodes, "/ResultsData/NodalDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, KratosHDF5.HDF5ConditionDataValueIO, self.model_part.GetCommunicator().LocalMesh().Conditions, "/ResultsData/ConditionDataValues/")
        self.__TestWriteVariableAndReadTensorAdaptor(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, KratosHDF5.HDF5ElementDataValueIO, self.model_part.GetCommunicator().LocalMesh().Elements, "/ResultsData/ElementDataValues/")

    @staticmethod
    def GetInputMDPAPath() -> Path:
        script_directory      = Path(__file__).absolute().parent
        kratos_root_directory = script_directory.parent.parent.parent
        test_input_directory  = kratos_root_directory / "kratos" / "tests" / "auxiliar_files_for_python_unittest"
        test_file_stem        = test_input_directory / "mdpa_files" / "test_processes"
        test_file_path        = Path(str(test_file_stem) + ".mdpa")

        if not test_file_path.is_file():
            raise FileNotFoundError("Test file not found: {}".format(test_file_path))

        return test_file_stem

if __name__ == "__main__":
    UnitTest.main()