from __future__ import print_function, absolute_import, division

import os
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MetisApplication as MetisApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestTrilinosRedistance(KratosUnittest.TestCase):

    def _ExpectedDistance(self,x,y,z):
        d = x
        if( d > 0.2):
            d = 0.2
        if( d < -0.2):
            d = -0.2
        return x
        #return -(math.sqrt(x**2+y**2+z**2) - 0.4)

    def tearDown(self):
        my_pid = self.model_part.GetCommunicator().MyPID()

        try:
            os.remove('coarse_sphere_' + my_pid + '.time')
            os.remove('coarse_sphere_' + my_pid + '.mdpa')
        except FileNotFoundError as e:
            pass

    def testTrilinosRedistance(self):
        # Set the model part
        self.model_part = KratosMultiphysics.ModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        # Import the model part, perform the partitioning and create communicators
        import_settings = KratosMultiphysics.Parameters("""{
            "echo_level" : 0,
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" : "coarse_sphere"
            }
        }""")
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.model_part, import_settings)
        TrilinosModelPartImporter.ExecutePartitioningAndReading()
        TrilinosModelPartImporter.CreateCommunicators()

        # Recall to set the buffer size
        self.model_part.SetBufferSize(2)

        # Initialize the DISTANCE values
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0, self._ExpectedDistance(node.X,node.Y,node.Z))

        # Set the output utility TODO: Remove after debugging
        import gid_output_process_mpi
        out_params = KratosMultiphysics.Parameters("""{
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostBinary",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteElementsOnly",
                    "MultiFileFlag": "SingleFile"
                },
                "file_label": "time",
                "output_control_type": "step",
                "output_frequency": 1.0,
                "body_output": true,
                "node_output": false,
                "skin_output": false,
                "plane_output": [],
                "nodal_results": ["DISTANCE"],
                "nodal_nonhistorical_results": [],
                "nodal_flags_results": [],
                "gauss_point_results": [],
                "additional_list_files": []
            },
            "point_data_configuration": []
        }""")
        gid_io = gid_output_process_mpi.GiDOutputProcessMPI(self.model_part, "coarse_sphere", out_params)
        gid_io.ExecuteInitialize()
        gid_io.ExecuteBeforeSolutionLoop()
        
        # Fake time advance
        self.model_part.CloneTimeStep(0.0)
        gid_io.ExecuteInitializeSolutionStep()
        gid_io.ExecuteFinalizeSolutionStep()
        if (gid_io.IsOutputStep()):
            gid_io.PrintOutput()
        
        self.model_part.CloneTimeStep(1.0)
        gid_io.ExecuteInitializeSolutionStep()

        # Set the utility and compute the variational distance values
        import trilinos_linear_solver_factory
        trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "AmesosSolver" }"""))

        epetra_comm = TrilinosApplication.CreateCommunicator()

        max_iterations = 2

        TrilinosApplication.TrilinosVariationalDistanceCalculationProcess3D(
            epetra_comm, self.model_part, trilinos_linear_solver, max_iterations).Execute()

        # Print results TODO: Remove after debugging
        gid_io.ExecuteFinalizeSolutionStep()
        if (gid_io.IsOutputStep()):
            gid_io.PrintOutput()

        self.model_part.CloneTimeStep(2.0)
        gid_io.ExecuteInitializeSolutionStep()
        gid_io.ExecuteFinalizeSolutionStep()
        if (gid_io.IsOutputStep()):
            gid_io.PrintOutput()
        gid_io.ExecuteFinalize()

        self.model_part.GetCommunicator().Barrier()

        # Check the obtained values
        max_distance = -1.0;
        min_distance = +1.0
        for node in self.model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        self.model_part.GetCommunicator().MinAll(min_distance)
        self.model_part.GetCommunicator().MaxAll(max_distance)

        self.assertAlmostEqual(max_distance, 0.44556526310761013)
        self.assertAlmostEqual(min_distance,-0.504972246827639)
        
        
if __name__ == '__main__':
    test = TestTrilinosRedistance()
    test.testTrilinosRedistance()
