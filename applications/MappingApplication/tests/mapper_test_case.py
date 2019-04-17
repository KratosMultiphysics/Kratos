from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
from KratosMultiphysics import KratosUnittest
data_comm = KM.DataCommunicator.GetDefault()
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "mdpa_files", fileName)

class MapperTestCase(KratosUnittest.TestCase):
    """This test case is designed for performing multiple test with the same modelparts
    this way the partitioning has to be done only once
    The values in the ModelParts are re-initialized after every
    """
    @classmethod
    def setUpModelParts(cls, mdpa_file_name_origin, mdpa_file_name_destination):
        cls.current_model = KM.Model()
        cls.model_part_origin = cls.current_model.CreateModelPart("origin")
        cls.model_part_destination = cls.current_model.CreateModelPart("destination")

        cls.input_file_origin      = GetFilePath(mdpa_file_name_origin)
        cls.input_file_destination = GetFilePath(mdpa_file_name_destination)

        cls.current_model = KM.Model()
        cls.model_part_origin = cls.current_model.CreateModelPart("origin")
        cls.model_part_destination = cls.current_model.CreateModelPart("destination")

        # list of variables involved in the Mapper-Tests
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.PRESSURE)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.FORCE)

        cls.model_part_destination.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.VELOCITY)

        cls.model_part_origin.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_destination.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_origin.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_origin.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes

        if data_comm.IsDistributed():
            ReadDistributedModelPart(cls.model_part_origin, cls.input_file_origin)
            ReadDistributedModelPart(cls.model_part_destination, cls.input_file_destination)
        else:
            ReadModelPart(cls.model_part_origin, cls.input_file_origin)
            ReadModelPart(cls.model_part_destination, cls.input_file_destination)

    def setUp(self):
        # reset the ModelPart
        # initialize it with random values, such that I am not accidentially
        # checking against 0.0
        default_scalar = -12345.6789
        default_vector = KM.Vector([111.222, -222.999, 333.444])

        for node in self.model_part_origin.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, default_scalar)
            node.SetSolutionStepValue(KM.FORCE, default_vector)
            node.SetValue(KM.PRESSURE, default_scalar)
            node.SetValue(KM.FORCE, default_vector)

        for element in self.model_part_origin.Elements:
            element.SetValue(KM.PRESSURE, default_scalar)
            element.SetValue(KM.FORCE, default_vector)

        for condition in self.model_part_origin.Conditions:
            condition.SetValue(KM.PRESSURE, default_scalar)
            condition.SetValue(KM.FORCE, default_vector)


        for node in self.model_part_destination.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, default_scalar)
            node.SetSolutionStepValue(KM.VELOCITY, default_vector)
            node.SetValue(KM.TEMPERATURE, default_scalar)
            node.SetValue(KM.VELOCITY, default_vector)

        for element in self.model_part_destination.Elements:
            element.SetValue(KM.TEMPERATURE, default_scalar)
            element.SetValue(KM.VELOCITY, default_vector)

        for condition in self.model_part_destination.Conditions:
            condition.SetValue(KM.TEMPERATURE, default_scalar)
            condition.SetValue(KM.VELOCITY, default_vector)

def ReadModelPart(model_part, mdpa_file_name):
    import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER

    KM.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

def ReadDistributedModelPart(model_part, mdpa_file_name):
    from KratosMultiphysics.TrilinosApplication import trilinos_import_model_part_utility
    model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

    importer_settings = KM.Parameters("""{
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": \"""" + mdpa_file_name + """\",
            "partition_in_memory" : true
        },
        "echo_level" : 0
    }""")

    model_part_import_util = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(model_part, importer_settings)
    model_part_import_util.ImportModelPart()
    model_part_import_util.CreateCommunicators()
