import pathlib

from KratosMultiphysics import KratosUnittest, Model, ModelPartIO, PARTITION_INDEX, FRACTIONAL_STEP, Logger
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.MetisApplication import MetisDivideHeterogeneousInputProcess

class QuadPartitionTest(KratosUnittest.TestCase):

    def testConditionPartitioning(self):
        test_file = "quads"
        partitions = 3
        sync_conditions = True
        domain_size = 2
        verbosity = 0

        work_dir = pathlib.Path(__file__).parent.absolute()
        model = Model()
        io_flags = ModelPartIO.READ|ModelPartIO.SKIP_TIMER|ModelPartIO.IGNORE_VARIABLES_ERROR

        # Register clean up operation for files created by the partitioner
        self.addCleanup(
            DeleteDirectoryIfExisting,
            work_dir / f'{test_file}_partitioned')

        with KratosUnittest.WorkFolderScope(work_dir, ''):
            io = ModelPartIO(test_file, io_flags)

            partitioner = MetisDivideHeterogeneousInputProcess(
                io, partitions, domain_size, verbosity, sync_conditions)
            partitioner.Execute()

            for i in range(partitions):
                model_part = model.CreateModelPart(f'test_{i}')
                model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

                file_path = f'{test_file}_partitioned/{test_file}_{i}'
                partitioned_io = ModelPartIO(file_path, io_flags)
                partitioned_io.ReadModelPart(model_part)

                for condition in model_part.Conditions:
                    parent = condition.GetValue(FRACTIONAL_STEP)
                    self.assertTrue(parent in model_part.Elements)

if __name__ == '__main__':
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
