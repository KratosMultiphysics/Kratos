import KratosMultiphysics as Kratos
from KratosMultiphysics import KratosUnittest

class SerialAndDistributedTestCase(KratosUnittest.TestCase):
    """This test case is designed for performing multiple test with the same modelparts,
    this way the partitioning has to be done only once
    The values in the ModelParts need to be re-initialized after every test (e.g. in setUp)
    """
    @classmethod
    def setUpClass(self):
        cls.model = Kratos.Model()

    @classmethod
    def _ReadModelPart(cls, mdpa_file_name, model_part):
        """This method is designed to read a ModelPart.

        The ModelPart is filled with respective nodes, conditions, elements
        from mdpa file using either MPI or Serial reading depending on run type.

        Args:
            mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
            model_part (Kratos.ModelPart): ModelPart to be filled
        """
        if Kratos.DOMAIN_SIZE not in model_part.ProcessInfo:
            raise Exception('"PROCESS_INFO" needs to be specified!')

        if model_part.NumberOfNodes() > 0:
            raise Exception("ModelPart must no contain Nodes!")

        communicator = Kratos.DataCommunicator.GetDefault()
        if communicator.IsDistributed():
            cls.__ReadDistributedModelPart(mdpa_file_name, model_part)
        else:
            cls.__ReadModelPart(mdpa_file_name, model_part)

    @classmethod
    def __ReadModelPart(cls, mdpa_file_name, model_part):
        """Reads mdpa file

        This method reads mdpa file and fills given model_part accordingly without MPI

        Args:
            mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
            model_part (Kratos.ModelPart): ModelPart to be filled
        """
        import_flags = Kratos.ModelPartIO.READ | Kratos.ModelPartIO.SKIP_TIMER
        Kratos.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

    @classmethod
    def __ReadDistributedModelPart(cls, mdpa_file_name, model_part):
        """Reads mdpa file

        This method reads mdpa file and fills given model_part accordingly using MPI

        Args:
            mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
            model_part (Kratos.ModelPart): ModelPart to be filled
        """
        KratosUnittest.skipIfApplicationsNotAvailable("MetisApplication")

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
