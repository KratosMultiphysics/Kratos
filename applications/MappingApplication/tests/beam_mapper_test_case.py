import KratosMultiphysics as KM
from KratosMultiphysics import KratosUnittest
import os
from KratosMultiphysics.testing import utilities as testing_utils

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "mdpa_files", fileName)

class BeamMapperTestCase(KratosUnittest.TestCase):
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
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.MOMENT)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.ROTATION)

        cls.model_part_destination.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.VELOCITY)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.REACTION)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT)

        cls.model_part_origin.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_destination.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_origin.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_origin.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes

        testing_utils.ReadModelPart(cls.input_file_origin, cls.model_part_origin)
        testing_utils.ReadModelPart(cls.input_file_destination, cls.model_part_destination)

    def setUp(self):
        # reset the ModelPart
        # initialize it with random values, such that I am not accidentially
        # checking against 0.0
        default_scalar = -12345.6789
        default_vector = KM.Vector([111.222, -222.999, 333.444])

        for node in self.model_part_origin.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, default_scalar)
            node.SetSolutionStepValue(KM.FORCE, default_vector)
            node.SetSolutionStepValue(KM.MOMENT, default_vector)
            node.SetValue(KM.PRESSURE, default_scalar)
            node.SetValue(KM.FORCE, default_vector)
            node.SetValue(KM.MOMENT, default_vector)

        for element in self.model_part_origin.Elements:
            element.SetValue(KM.PRESSURE, default_scalar)
            element.SetValue(KM.FORCE, default_vector)
            element.SetValue(KM.MOMENT, default_vector)

        for condition in self.model_part_origin.Conditions:
            condition.SetValue(KM.PRESSURE, default_scalar)
            condition.SetValue(KM.FORCE, default_vector)
            condition.SetValue(KM.MOMENT, default_vector)

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