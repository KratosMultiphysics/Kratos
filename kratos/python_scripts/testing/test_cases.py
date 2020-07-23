import KratosMultiphysics as Kratos
from KratosMultiphysics import KratosUnittest

class KratosCoreTestCase(KratosUnittest.TestCase):
    """This test case is designed for performing multiple test with the same modelparts,
    this way the partitioning has to be done only once
    The values in the ModelParts are re-initialized after every test
    """
    @classmethod
    def setUpClass(self):
        cls.model = Kratos.Model()

    @classmethod
    def _ReadModelPart(cls, mdpa_file_name, model_part):
        """This method setup the model part once for all the tests in class.

        This method read creates a modelpart by the given name with given domain size,
        buffer size. Then created model part is filled with respective nodes, conditions, elements
        from mdpa file using either MPI or Serial reading depending on run type.

        Args:
            model_part: Name of test model part
        """
        if not Kratos.DOMAIN_SIZE in model_part.ProcessInfo:
            raise Exception('"PROCESS_INFO" needs to be specified!')

        if model_part.NumberOfNodes() > 0:
            raise Exception("ModelPart must no contain Nodes!")

        communicator = Kratos.DataCommunicator.GetDefault()
        if communicator.IsDistributed():
            cls.__ReadDistributedModelPart(mdpa_file_name)
        else:
            cls.__ReadModelPart(mdpa_file_name)

    @classmethod
    def __ReadModelPart(cls, mdpa_file_name, model_part):
        """Reads mdpa file

        This method reads mdpa file and fills given model_part accordingly without MPI

        Args:
            model_part (Kratos.ModelPart): Model part to be filled
            mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
        """
        import_flags = Kratos.ModelPartIO.READ | Kratos.ModelPartIO.SKIP_TIMER
        Kratos.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

    @classmethod
    def __ReadDistributedModelPart(cls, mdpa_file_name, model_part):
        """Reads mdpa file

        This method reads mdpa file and fills given model_part accordingly using MPI

        Args:
            model_part (Kratos.ModelPart): Model part to be filled
            mdpa_file_name (str): Name of the mdpa file (without '.mdpa' extension)
        """
        from KratosMultiphysics.mpi import distributed_import_model_part_utility
        model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)

        importer_settings = Kratos.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + mdpa_file_name + """\",
                "partition_in_memory" : true
            },
            "echo_level" : 0
        }""")

        model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, importer_settings)
        model_part_import_util.ImportModelPart()
        model_part_import_util.CreateCommunicators()
