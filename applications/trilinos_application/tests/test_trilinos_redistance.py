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

        # Remove the .time file
        if (my_pid == 0):
            try:
                os.remove('coarse_sphere.time')
            except FileNotFoundError as e:
                pass

        # Remove the Metis partitioning files
        try:
            os.remove('coarse_sphere_' + str(my_pid) + '.time')
            os.remove('coarse_sphere_' + str(my_pid) + '.mdpa')
        except FileNotFoundError as e:
            pass

    def testTrilinosRedistance(self):
        # Set the model part
        current_model = KratosMultiphysics.Model()
        self.model_part = current_model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
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
        TrilinosModelPartImporter.ImportModelPart()
        TrilinosModelPartImporter.CreateCommunicators()

        # Recall to set the buffer size
        self.model_part.SetBufferSize(2)

        # Initialize the DISTANCE values
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0, self._ExpectedDistance(node.X,node.Y,node.Z))

        # Fake time advance
        self.model_part.CloneTimeStep(1.0)

        # Set the utility and compute the variational distance values
        import trilinos_linear_solver_factory
        trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "AmesosSolver" }"""))

        epetra_comm = TrilinosApplication.CreateCommunicator()

        max_iterations = 2
        TrilinosApplication.TrilinosVariationalDistanceCalculationProcess3D(
            epetra_comm, self.model_part, trilinos_linear_solver, max_iterations).Execute()

        # Check the obtained values
        max_distance = -1.0
        min_distance = +1.0
        for node in self.model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        min_distance = self.model_part.GetCommunicator().MinAll(min_distance)
        max_distance = self.model_part.GetCommunicator().MaxAll(max_distance)

        self.assertAlmostEqual(max_distance, 0.44556526310761013) # Serial max_distance
        self.assertAlmostEqual(min_distance,-0.504972246827639) # Serial min_distance

if __name__ == '__main__':
    KratosUnittest.main()
