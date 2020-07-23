import KratosMultiphysics as Kratos
from KratosMultiphysics import KratosUnittest

class KratosCoreTestCase(KratosUnittest.TestCase):
    """This test case is designed for performing multiple test with the same modelparts
    this way the partitioning has to be done only once
    The values in the ModelParts are re-initialized after every test
    """
    @classmethod
    def setUpModelParts(cls, model_part_name, mdpa_file_name, domain_size = 2, buffer_size = 1):
        """This method setup the model part once for all the tests in class.

        This method read creates a modelpart by the given name with given domain size,
        buffer size. Then created model part is filled with respective nodes, conditions, elements
        from mdpa file using either MPI or Serial reading depending on run type.

        Args:
            model_part_name (str): Name of test model part
            mdpa_file_name (str): Name of mdpa file name (without '.mdpa' extension)
            domain_size (int): Domain size of model part
            buffer_size (int): Buffer size required by model part
        """
        cls.current_model = Kratos.Model()
        cls.model_part = cls.current_model.CreateModelPart(model_part_name)
        cls.model_part.SetBufferSize(buffer_size)
        cls.model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, domain_size)

        cls.AddVariables()

        communicator = Kratos.DataCommunicator.GetDefault()
        if communicator.IsDistributed():
            cls.__ReadDistributedModelPart(mdpa_file_name)
        else:
            cls.__ReadModelPart(mdpa_file_name)

    @classmethod
    def AddVariables(cls):
        """Add variables to solution step data value container
        """
        pass

    @classmethod
    def GetModelPart(cls):
        """Returns created model part

        Returns:
            Kratos.ModelPart: Returns created model part
        """
        return cls.model_part

    @classmethod
    def GetModel(cls):
        """Returns created model part

        Returns:
            Kratos.Model: Returns created model
        """
        return cls.current_model

    @classmethod
    def __ReadModelPart(cls, mdpa_file_name):
        """Reads mdpa file

        This method reads mdpa file and fills given model_part accordingly without MPI

        Args:
            model_part (Kratos.ModelPart): Model part to be filled
            mdpa_file_name (str): Name of the mdpa file (without ".mdpa" extension)
        """
        import_flags = Kratos.ModelPartIO.READ | Kratos.ModelPartIO.SKIP_TIMER
        Kratos.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(cls.GetModelPart())

    @classmethod
    def __ReadDistributedModelPart(cls, mdpa_file_name):
        """Reads mdpa file

        This method reads mdpa file and fills given model_part accordingly using MPI

        Args:
            model_part (Kratos.ModelPart): Model part to be filled
            mdpa_file_name (str): Name of the mdpa file (without '.mdpa' extension)
        """
        from KratosMultiphysics.mpi import distributed_import_model_part_utility
        cls.GetModelPart().AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)

        importer_settings = Kratos.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + mdpa_file_name + """\",
                "partition_in_memory" : true
            },
            "echo_level" : 0
        }""")

        model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(cls.GetModelPart(), importer_settings)
        model_part_import_util.ImportModelPart()
        model_part_import_util.CreateCommunicators()