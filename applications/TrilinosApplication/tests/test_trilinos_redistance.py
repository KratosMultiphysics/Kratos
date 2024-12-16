import pathlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.kratos_utilities as KratosUtils

from KratosMultiphysics.testing.utilities import ReadDistributedModelPart
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

def GetFilePath(fileName):
    return str(pathlib.Path(__file__).absolute().parent / fileName)

class TestTrilinosRedistance(KratosUnittest.TestCase):

    @classmethod
    def _ExpectedDistance(self,x,y,z):
        d = x
        if( d > 0.2):
            d = 0.2
        if( d < -0.2):
            d = -0.2
        return x
        #return -(math.sqrt(x**2+y**2+z**2) - 0.4)

    @classmethod
    def _ExpectedLinearDistance(self,x,x_0):
        return x - x_0

    def setUp(self):
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("Main",2)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ReadDistributedModelPart(GetFilePath( "auxiliary_files/mdpa_files/coarse_sphere"), self.model_part)

    def testTrilinosRedistance(self):
        # Initialize the DISTANCE values
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0, self._ExpectedDistance(node.X,node.Y,node.Z))

        # Fake time advance
        self.model_part.CloneTimeStep(1.0)

        # Set the utility and compute the variational distance values
        trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "amesos" }"""))

        epetra_comm = TrilinosApplication.CreateEpetraCommunicator(self.model_part.GetCommunicator().GetDataCommunicator())

        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "Main",
            "max_iterations" : 2
        }""")
        TrilinosApplication.TrilinosVariationalDistanceCalculationProcess3D(epetra_comm, self.current_model, trilinos_linear_solver, settings).Execute()

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

        self.assertAlmostEqual(max_distance, 0.44556526310761013) # Serial max_distance
        self.assertAlmostEqual(min_distance,-0.504972246827639) # Serial min_distance

    def testTrilinosParallelRedistance(self):
        # Initialize the DISTANCE values
        x_zero_dist = 0.0
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0, self._ExpectedLinearDistance(node.X, x_zero_dist))

        # Fake time advance
        self.model_part.CloneTimeStep(1.0)

        # Calculate NODAL_AREA
        domain_size = 3
        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(self.model_part, domain_size)
        nodal_area_process.Execute()

        # Set the parallel distance calculator
        settings = KratosMultiphysics.Parameters("""{
            "max_levels" : 10,
            "max_distance" : 100.0,
            "calculate_exact_distances_to_plane" : true
        }""")
        distance_calculator = KratosMultiphysics.ParallelDistanceCalculationProcess3D(
            self.model_part,
            settings)
        distance_calculator.Execute()

        # Check the obtained values
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), self._ExpectedLinearDistance(node.X, x_zero_dist), 10)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
