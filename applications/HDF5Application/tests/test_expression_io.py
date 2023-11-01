import numpy as np
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

class TestExpressionIO(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        ReadModelPart(str(TestExpressionIO.GetInputMDPAPath()), cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, node.Id)
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))
            node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, Kratos.Vector(4, node.Id))
            node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix(2, 3, node.Id))

        for condition in cls.model_part.Conditions:
            condition.SetValue(Kratos.PRESSURE, condition.Id)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id, condition.Id + 1, condition.Id + 2]))
            condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, Kratos.Vector(4, condition.Id))
            condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix(2, 3, condition.Id))

        for element in cls.model_part.Elements:
            element.SetValue(Kratos.PRESSURE, element.Id)
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([element.Id, element.Id + 1, element.Id + 2]))
            element.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, Kratos.Vector(4, element.Id))
            element.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix(2, 3, element.Id))

        cls.data_communicator: Kratos.DataCommunicator = cls.model_part.GetCommunicator().GetDataCommunicator()

    def setUp(self) -> None:
        if self.data_communicator.IsDistributed():
            parameters = Kratos.Parameters("""{
                "file_name" : "test_expression_io.h5",
                "file_access_mode": "truncate",
                "file_driver": "mpio",
                "echo_level" : 0
            }""")
        else:
            parameters = Kratos.Parameters("""{
                "file_name" : "test_expression_io.h5",
                "file_access_mode": "truncate",
                "file_driver": "core",
                "echo_level" : 0
            }""")
        self.h5_file = KratosHDF5.HDF5File(self.data_communicator, parameters)

    def tearDown(self) -> None:
        self.data_communicator.Barrier()
        DeleteFileIfExisting("test_expression_io.h5")

    def __TestContainerExpressions(self, variable):
        attribs_in = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")
        nodal_expression = Kratos.Expression.NodalExpression(self.model_part)
        condition_expression = Kratos.Expression.ConditionExpression(self.model_part)
        element_expression = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(nodal_expression, variable, False)
        Kratos.Expression.VariableExpressionIO.Read(condition_expression, variable)
        Kratos.Expression.VariableExpressionIO.Read(element_expression, variable)
        expression_io = KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix": "/expressions"}"""), self.h5_file)
        expression_io.Write("nodal", nodal_expression, attribs_in)
        expression_io.Write("condition", condition_expression, attribs_in)
        expression_io.Write("element", element_expression, attribs_in)

        nodal_expression_out = Kratos.Expression.NodalExpression(self.model_part)
        condition_expression_out = Kratos.Expression.ConditionExpression(self.model_part)
        element_expression_out = Kratos.Expression.ElementExpression(self.model_part)
        nodal_attribs_out = expression_io.Read("nodal", nodal_expression_out)
        condition_attribs_out = expression_io.Read("condition", condition_expression_out)
        element_attribs_out = expression_io.Read("element", element_expression_out)
        self.assertTrue(nodal_attribs_out.IsEquivalentTo(attribs_in))
        self.assertTrue(condition_attribs_out.IsEquivalentTo(attribs_in))
        self.assertTrue(element_attribs_out.IsEquivalentTo(attribs_in))
        self.assertEqual(np.linalg.norm((nodal_expression_out - nodal_expression).Evaluate()), 0)
        self.assertEqual(np.linalg.norm((condition_expression_out - condition_expression).Evaluate()), 0)
        self.assertEqual(np.linalg.norm((element_expression_out - element_expression).Evaluate()), 0)

    def test_ReadWriteExpression(self):
        pass
        # TODO: this test needs to be enabled when then the problem of return type deduction in python is solved which is caused by the use of Trampoline class for expressions.
        # nodal_expression = Kratos.Expression.NodalExpression(self.model_part)
        # Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.PRESSURE, False)
        # expression_io = KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix": "/expressions"}"""), self.h5_file)

        # in_attribs = Kratos.Parameters("""{"custom_attrib": "custom_value"}""")
        # expression_io.Write("test", nodal_expression.GetExpression(), in_attribs)
        # out_expression, out_attribs = expression_io.Read("test")

    def __TestWriteVariableAndReadExpression(self, variable):
        nodal_params = Kratos.Parameters("""{"prefix": "/ResultsData" }""")
        nodal_params.AddStringArray("list_of_variables", [variable.Name()])
        nodal_io = KratosHDF5.HDF5NodalDataValueIO(nodal_params, self.h5_file)
        nodal_io.Write(self.model_part)
        nodal_expression_orig = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(nodal_expression_orig, variable, False)
        nodal_expression_read = Kratos.Expression.NodalExpression(self.model_part)
        KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix": "/ResultsData/NodalDataValues/"}"""), self.h5_file).Read(variable.Name(), nodal_expression_read)
        self.assertEqual(np.linalg.norm((nodal_expression_read - nodal_expression_orig).Evaluate()), 0)

        condition_params = Kratos.Parameters("""{"prefix": "/ResultsData" }""")
        condition_params.AddStringArray("list_of_variables", [variable.Name()])
        condition_io = KratosHDF5.HDF5ConditionDataValueIO(condition_params, self.h5_file)
        condition_io.Write(self.model_part)
        condition_expression_orig = Kratos.Expression.ConditionExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(condition_expression_orig, variable)
        condition_expression_read = Kratos.Expression.ConditionExpression(self.model_part)
        KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix": "/ResultsData/ConditionDataValues/"}"""), self.h5_file).Read(variable.Name(), condition_expression_read)
        self.assertEqual(np.linalg.norm((condition_expression_read - condition_expression_orig).Evaluate()), 0)

        element_params = Kratos.Parameters("""{"prefix": "/ResultsData" }""")
        element_params.AddStringArray("list_of_variables", [variable.Name()])
        element_io = KratosHDF5.HDF5ElementDataValueIO(element_params, self.h5_file)
        element_io.Write(self.model_part)
        element_expression_orig = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(element_expression_orig, variable)
        element_expression_read = Kratos.Expression.ElementExpression(self.model_part)
        KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix": "/ResultsData/ElementDataValues/"}"""), self.h5_file).Read(variable.Name(), element_expression_read)
        self.assertEqual(np.linalg.norm((element_expression_read - element_expression_orig).Evaluate()), 0)

    def test_ReadWriteExpressionScalar(self):
        self.__TestContainerExpressions(Kratos.PRESSURE)

    def test_ReadWriteExpressionArray3(self):
        self.__TestContainerExpressions(Kratos.VELOCITY)

    def test_ReadWriteExpressionVector(self):
        self.__TestContainerExpressions(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR)

    def test_ReadWriteExpressionMatrix(self):
        self.__TestContainerExpressions(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

    def test_WriteVariableReadExpressionScalar(self):
        self.__TestWriteVariableAndReadExpression(Kratos.PRESSURE)

    def test_WriteVariableReadExpressionArray3(self):
        self.__TestWriteVariableAndReadExpression(Kratos.VELOCITY)

    def test_WriteVariableReadExpressionVector(self):
        self.__TestWriteVariableAndReadExpression(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR)

    def test_WriteVariableReadExpressionMatrix(self):
        self.__TestWriteVariableAndReadExpression(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)


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