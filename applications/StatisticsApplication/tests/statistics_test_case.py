from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
from KratosMultiphysics import KratosUnittest

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

        communicator = Kratos.DataCommunicator.GetDefault()
        if communicator.IsDistributed():
            ReadDistributedModelPart(cls.model_part, mdpa_file_name)
        else:
            ReadModelPart(cls.model_part, mdpa_file_name)

    @classmethod
    def AddVariables(cls):
        pass

    def setUp(self):
        InitializeModelPartVariables(self.model_part)

    def GetModelPart(self):
        return self.model_part

    def GetModel(self):
        return self.current_model

def ReadModelPart(model_part, mdpa_file_name):
    import_flags = Kratos.ModelPartIO.READ | Kratos.ModelPartIO.SKIP_TIMER
    Kratos.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

def ReadDistributedModelPart(model_part, mdpa_file_name):
    from KratosMultiphysics.mpi import distributed_import_model_part_utility
    model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)

    importer_settings = Kratos.Parameters("""{
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": \"""" + mdpa_file_name + """\",
            "partition_in_memory" : true
        },
        "echo_level" : 0
    }""")

    model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, importer_settings)
    model_part_import_util.ImportModelPart()
    model_part_import_util.CreateCommunicators()