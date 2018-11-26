from __future__ import print_function, absolute_import, division

import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MetisApplication as MetisApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.kratos_utilities as KratosUtils

import trilinos_import_model_part_utility
import trilinos_linear_solver_factory

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseDistance(x, y, z):
    if (x <= 5.0):
        return -0.16*x**2 + 0.8*x
    else:
        return 0.0

def ConvectionVelocity(x, y, z):
    vel = KratosMultiphysics.Vector(3, 0.0)
    vel[0] = 1.0
    return vel

class TestTrilinosLevelSetConvection(KratosUnittest.TestCase):

    def setUp(self):
        self.parameters = """{
            "echo_level" : 0,
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" : "levelset_convection_process_mesh"
            }
        } """

    def tearDown(self):
        my_pid = self.model_part.GetCommunicator().MyPID()

        # Remove the .time file
        if (my_pid == 0):
            KratosUtils.DeleteFileIfExisting("levelset_convection_process_mesh.time")

        # Remove the Metis partitioning files
        KratosUtils.DeleteFileIfExisting("levelset_convection_process_mesh_" + str(my_pid) + ".time")
        KratosUtils.DeleteFileIfExisting("levelset_convection_process_mesh_" + str(my_pid) + ".mdpa")

        # While compining in debug, in memory partitioner also writes down the mpda in plain text
        # and needs to be cleaned.
        KratosUtils.DeleteFileIfExisting("debug_modelpart_" + str(my_pid) + ".mdpa")

    def test_trilinos_levelset_convection(self):
        current_model = KratosMultiphysics.Model()
        self.model_part = current_model.CreateModelPart("Main",2)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)


        # Import the model part, perform the partitioning and create communicators
        import_settings = KratosMultiphysics.Parameters(self.parameters)
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.model_part, import_settings)
        TrilinosModelPartImporter.ImportModelPart()
        TrilinosModelPartImporter.CreateCommunicators()

        # Recall to set the buffer size
        self.model_part.SetBufferSize(2)

        # Set the initial distance field and the convection velocity
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, BaseDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, ConvectionVelocity(node.X,node.Y,node.Z))

        # Fix the left side values
        for node in self.model_part.Nodes:
            if node.X < 0.001:
                node.Fix(KratosMultiphysics.DISTANCE)

        # Set the Trilinos linear solver and Epetra communicator
        trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "AmesosSolver" }""")
        )

        epetra_comm = TrilinosApplication.CreateCommunicator()

        # Fake time advance
        self.model_part.CloneTimeStep(40.0)

        # Convect the distance field
        TrilinosApplication.TrilinosLevelSetConvectionProcess2D(
            epetra_comm,
            KratosMultiphysics.DISTANCE,
            self.model_part,
            trilinos_linear_solver).Execute()

        # Check the obtained values
        max_distance = -1.0
        min_distance = +1.0

        for node in self.model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        min_distance = self.model_part.GetCommunicator().MinAll(min_distance)
        max_distance = self.model_part.GetCommunicator().MaxAll(max_distance)

        self.assertAlmostEqual(max_distance, 0.7904255118014996)
        self.assertAlmostEqual(min_distance,-0.11710292469993094)

class TestTrilinosLevelSetConvectionInMemory(TestTrilinosLevelSetConvection):

    def setUp(self):
        self.parameters = """{
            "echo_level" : 0,
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" : "levelset_convection_process_mesh",
                "partition_in_memory" : true
            }
        } """

if __name__ == '__main__':
    KratosUnittest.main()