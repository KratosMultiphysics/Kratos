from __future__ import print_function, absolute_import, division

import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MetisApplication as MetisApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.kratos_utilities as KratosUtils

from KratosMultiphysics.mpi import distributed_import_model_part_utility
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

from KratosMultiphysics import ParallelEnvironment

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseDistance(x, y, z):
    if (x <= 5.0):
        return -0.16*x**2 + 0.8*x
    else:
        return 0.0

def BaseJumpedDistance(x, y, z):
    if (x >= 5.0 and x <= 15.0):
        return 1.0
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
                "input_filename" : \"""" + GetFilePath("levelset_convection_process_mesh") + """\"
            }
        } """

    def tearDown(self):
        my_pid = self.model_part.GetCommunicator().MyPID()

        # Remove the .time file
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
        DistributedModelPartImporter = distributed_import_model_part_utility.DistributedImportModelPartUtility(self.model_part, import_settings)
        DistributedModelPartImporter.ImportModelPart()
        DistributedModelPartImporter.CreateCommunicators()

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
            KratosMultiphysics.Parameters("""{"solver_type" : "amesos" }""")
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

        comm = self.model_part.GetCommunicator().GetDataCommunicator()
        min_distance = comm.MinAll(min_distance)
        max_distance = comm.MaxAll(max_distance)

        self.assertAlmostEqual(max_distance, 0.7333041045431626)
        self.assertAlmostEqual(min_distance,-0.06371359024393104)

    def test_trilinos_levelset_convection_BFECC(self):
        current_model = KratosMultiphysics.Model()
        self.model_part = current_model.CreateModelPart("Main",2)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        # Import the model part, perform the partitioning and create communicators
        import_settings = KratosMultiphysics.Parameters(self.parameters)
        DistributedModelPartImporter = distributed_import_model_part_utility.DistributedImportModelPartUtility(self.model_part, import_settings)
        DistributedModelPartImporter.ImportModelPart()
        DistributedModelPartImporter.CreateCommunicators()

        # Recall to set the buffer size
        self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # Set the initial distance field and the convection velocity
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, BaseJumpedDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, ConvectionVelocity(node.X,node.Y,node.Z))

        # Fix the left side values
        for node in self.model_part.Nodes:
            if node.X < 0.001:
                node.Fix(KratosMultiphysics.DISTANCE)

        # Set the Trilinos linear solver and Epetra communicator
        trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "amesos" }""")
        )
        epetra_comm = TrilinosApplication.CreateCommunicator()
        comm = ParallelEnvironment.GetDefaultDataCommunicator()
        #self.model_part.GetCommunicator().GetDataCommunicator()

        # Fake time advance
        self.model_part.CloneTimeStep(30.0)

        #kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(
                comm, self.model_part).Execute()

        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
            self.model_part,
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "levelset_variable_name" : "DISTANCE",
            "levelset_convection_variable_name" : "VELOCITY",
            "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
            "max_CFL" : 1.0,
            "max_substeps" : 0,
            "levelset_splitting" : false,
            "eulerian_error_compensation" : true,
            "cross_wind_stabilization_factor" : 0.7
        }""")
        TrilinosApplication.TrilinosLevelSetConvectionProcess2D(
            epetra_comm,
            self.model_part,
            trilinos_linear_solver,
            levelset_convection_settings).Execute()

        max_distance = -1.0
        min_distance = +1.0

        for node in self.model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        min_distance = comm.MinAll(min_distance)
        max_distance = comm.MaxAll(max_distance)

        # gid_output = GiDOutputProcess(model_part,
        #                            "levelset_test_2D",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","VELOCITY"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

        self.assertAlmostEqual(max_distance, 1.0617777301844604)
        self.assertAlmostEqual(min_distance, -0.061745786561321375)

class TestTrilinosLevelSetConvectionInMemory(TestTrilinosLevelSetConvection):

    def setUp(self):
        self.parameters = """{
            "echo_level" : 0,
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" : \"""" + GetFilePath("levelset_convection_process_mesh") + """\",
                "partition_in_memory" : true
            }
        } """

if __name__ == '__main__':
    KratosUnittest.main()