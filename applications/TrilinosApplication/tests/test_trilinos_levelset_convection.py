import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.kratos_utilities as KratosUtils

from KratosMultiphysics.testing.utilities import ReadDistributedModelPart
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory


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

    @classmethod
    def GetPartitioningParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + GetFilePath("levelset_convection_process_mesh") + """\",
                "partition_in_memory" : false
            },
            "echo_level" : 0
        }""")

    def setUp(self):
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("Main",2)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ReadDistributedModelPart(GetFilePath("levelset_convection_process_mesh"), self.model_part, self.GetPartitioningParameters())

    def tearDown(self):
        # Remove the Metis partitioning files
        KratosUtils.DeleteDirectoryIfExisting("levelset_convection_process_mesh_partitioned")
        # next test can only start after all the processes arrived here, otherwise race conditions with deleting the files can occur
        self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

    def test_trilinos_levelset_convection(self):
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

        epetra_comm = TrilinosApplication.CreateEpetraCommunicator(KratosMultiphysics.DataCommunicator.GetDefault())

        # Fake time advance
        self.model_part.CloneTimeStep(40.0)

        # Convect the distance field
        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "max_CFL" : 1.0,
            "max_substeps" : 0,
            "eulerian_error_compensation" : false,
            "element_type" : "levelset_convection_supg"
        }""")
        TrilinosApplication.TrilinosLevelSetConvectionProcess2D(
            epetra_comm,
            self.model_part,
            trilinos_linear_solver,
            levelset_convection_settings).Execute()

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
        epetra_comm = TrilinosApplication.CreateEpetraCommunicator(KratosMultiphysics.DataCommunicator.GetDefault())

        # Fake time advance
        self.model_part.CloneTimeStep(30.0)

        KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.model_part).Execute()

        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "max_CFL" : 1.0,
            "max_substeps" : 0,
            "eulerian_error_compensation" : true,
            "element_type" : "levelset_convection_supg"
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

        kratos_comm = self.model_part.GetCommunicator().GetDataCommunicator()
        min_distance = kratos_comm.MinAll(min_distance)
        max_distance = kratos_comm.MaxAll(max_distance)

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

        self.assertAlmostEqual(max_distance, 1.0634680107706003)
        self.assertAlmostEqual(min_distance, -0.06361967738862996)

class TestTrilinosLevelSetConvectionInMemory(TestTrilinosLevelSetConvection):

    @classmethod
    def GetPartitioningParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + GetFilePath("levelset_convection_process_mesh") + """\",
                "partition_in_memory" : true
            },
            "echo_level" : 0
        }""")

if __name__ == '__main__':
    KratosUnittest.main()
