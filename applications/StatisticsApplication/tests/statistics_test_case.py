import KratosMultiphysics as Kratos
from KratosMultiphysics import KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeModelPartVariables

class StatisticsTestCase(KratosUnittest.TestCase):
    """This test case is designed for performing multiple test with the same modelparts
    this way the partitioning has to be done only once
    The values in the ModelParts are re-initialized after every test
    """
    @classmethod
    def setUpModelParts(cls, mdpa_file_name):
        cls.current_model = Kratos.Model()
        cls.model_part = cls.current_model.CreateModelPart("test_model_part")
        cls.model_part.SetBufferSize(1)
        cls.model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)

        cls.AddVariables()
        with KratosUnittest.WorkFolderScope(".", __file__):
            ReadModelPart(mdpa_file_name, cls.model_part)

    @classmethod
    def AddVariables(cls):
        pass

    def setUp(self):
        InitializeModelPartVariables(self.model_part)

    def GetModelPart(self):
        return self.model_part

    def GetModel(self):
        return self.current_model
