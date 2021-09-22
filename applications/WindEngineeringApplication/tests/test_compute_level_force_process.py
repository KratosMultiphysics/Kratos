# External imports
import numpy

# Kratos imports
import KratosMultiphysics
from KratosMultiphysics.WindEngineeringApplication.compute_level_force_process import ComputeLevelForceProcess
from KratosMultiphysics.WindEngineeringApplication.compute_level_force_process import Factory as ComputeLevelForceProcessFactory
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils

# STL imports
import pathlib


def GetFilePath(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName


class TestComputeLevelForceProcess(UnitTest.TestCase):

    def GenerateModel(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Fill a unit cube with nodes and assign reactions to them
        coordinates = numpy.linspace(0.0, 1.0, num=6)
        node_counter = 0
        for x in coordinates:
            for y in coordinates:
                for z in coordinates:
                    node_counter += 1
                    node = model_part.CreateNewNode(node_counter, x, y, z)
                    node.SetSolutionStepValue(KratosMultiphysics.REACTION, [x*(y+z), -z, y])

        return model, model_part

    def testComputeLevelForceProcess(self):
        parameters = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "model_part_name"   : "main",
                "bottom_point"      : [0.0, 0.0, 0.0],
                "top_point"         : [0.9, 0.0, 0.0],
                "number_of_slabs"   : 3,
                "output_name_stub"  : ""
            }
        }
        """)
        parameters["Parameters"]["output_name_stub"].SetString(str(GetFilePath("test_compute_level_force_process_")))
        model, model_part = self.GenerateModel()
        ComputeLevelForceProcessFactory(
            parameters,
            model).ExecuteFinalizeSolutionStep()

        force, torque = ComputeLevelForceProcess.ParseOutput(GetFilePath("test_compute_level_force_process_0.dat"))
        self.assertAlmostEqual(force[0], 7.2)
        self.assertAlmostEqual(force[1], -36.0)
        self.assertAlmostEqual(force[2], 36.0)
        self.assertAlmostEqual(torque[0], 52.8)
        self.assertAlmostEqual(torque[1], 0.84)
        self.assertAlmostEqual(torque[2], -8.04)

        force, torque = ComputeLevelForceProcess.ParseOutput(GetFilePath("test_compute_level_force_process_1.dat"))
        self.assertAlmostEqual(force[0], 14.4)
        self.assertAlmostEqual(force[1], -18.0)
        self.assertAlmostEqual(force[2], 18.0)
        self.assertAlmostEqual(torque[0], 26.4)
        self.assertAlmostEqual(torque[1], 1.68)
        self.assertAlmostEqual(torque[2], -16.08)

        force, torque = ComputeLevelForceProcess.ParseOutput(GetFilePath("test_compute_level_force_process_2.dat"))
        self.assertAlmostEqual(force[0], 50.4)
        self.assertAlmostEqual(force[1], -36.0)
        self.assertAlmostEqual(force[2], 36.0)
        self.assertAlmostEqual(torque[0], 52.8)
        self.assertAlmostEqual(torque[1], 5.88)
        self.assertAlmostEqual(torque[2], -56.28)


    @classmethod
    def tearDown(cls):
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_compute_level_force_process_0.dat"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_compute_level_force_process_1.dat"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_compute_level_force_process_2.dat"))


if __name__ == "__main__":
    UnitTest.main() 