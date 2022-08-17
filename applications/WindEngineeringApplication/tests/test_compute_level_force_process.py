# Kratos imports
import KratosMultiphysics
from KratosMultiphysics.WindEngineeringApplication.compute_level_force_process import ComputeLevelForceProcess
from KratosMultiphysics.WindEngineeringApplication.compute_level_force_process import Factory as ComputeLevelForceProcessFactory
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.WindEngineeringApplication.test_case import SuiteFlags, TestCase
import KratosMultiphysics.kratos_utilities as KratosUtils

# STL imports
import pathlib


def GetFilePath(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName


class TestComputeLevelForceProcess(TestCase):
    suite_flags = [SuiteFlags.ALL, SuiteFlags.MPI]

    def GenerateModel(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Fill a unit cube with nodes and assign reactions to them
        coordinates = [index/5.0 for index in range(6)]
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
                "bottom_point"      : [-0.1, 0.0, 0.0],
                "top_point"         : [0.91, 0.0, 0.0],
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
        self.assertVectorAlmostEqual(force, [7.2, -36.0, 36.0])
        self.assertVectorAlmostEqual(torque, [52.8, 0.84, -8.04])

        force, torque = ComputeLevelForceProcess.ParseOutput(GetFilePath("test_compute_level_force_process_1.dat"))
        self.assertVectorAlmostEqual(force, [14.4, -18.0, 18.0])
        self.assertVectorAlmostEqual(torque, [26.4, 1.68, -16.08])

        force, torque = ComputeLevelForceProcess.ParseOutput(GetFilePath("test_compute_level_force_process_2.dat"))
        self.assertVectorAlmostEqual(force, [50.4, -36.0, 36.0])
        self.assertVectorAlmostEqual(torque, [52.8, 5.88, -56.28])


    @classmethod
    def tearDown(cls):
        KratosMultiphysics.Testing.GetDefaultDataCommunicator().Barrier()
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_compute_level_force_process_0.dat"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_compute_level_force_process_1.dat"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("test_compute_level_force_process_2.dat"))


if __name__ == "__main__":
    UnitTest.main()