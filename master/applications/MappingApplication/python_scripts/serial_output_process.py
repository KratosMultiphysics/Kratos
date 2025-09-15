import KratosMultiphysics as KM
from KratosMultiphysics import MappingApplication # registering the mappers
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput
from KratosMultiphysics.process_factory import KratosProcessFactory

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SerialOutputProcess(model, settings["Parameters"])


class SerialOutputProcess(KM.OutputProcess):
    """This process is used in distributed simulations to do post-processing on one rank
    For this first the results are mapped to a serial ModelPart, with which then postprocessing is done
    This is not the most efficient approach by principle, but sometimes it is necessary to collect
    all the results on one rank
    """
    def __init__(self, model, settings):
        super().__init__()

        self.settings = settings

        default_settings = KM.Parameters('''{
            "main_model_part_name_origin"      : "UNSPECIFIED",
            "main_model_part_name_destination" : "UNSPECIFIED",
            "mdpa_file_name_destination"  : "",
            "historical_variables_destination" : [],
            "destination_rank" : 0,
            "mapper_settings" :  {},
            "mapping_settings" : [],
            "output_process_settings" : {}
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        if len(settings["mapper_settings"].keys()) == 0:
            raise Exception('no "mapper_settings" were specified!')

        if settings["mapping_settings"].size() == 0:
            raise Exception('no "mapping_settings" were specified!')

        mdpa_file_name_destination = settings["mdpa_file_name_destination"].GetString()

        model_part_origin = model.GetModelPart(settings["main_model_part_name_origin"].GetString())
        if model_part_origin.IsSubModelPart():
            raise Exception('Origin ModelPart cannot be a SubModelPart!')

        self.model_part_destination = model.CreateModelPart(settings["main_model_part_name_destination"].GetString())
        if self.model_part_destination.IsSubModelPart():
            raise Exception('Destination ModelPart cannot be a SubModelPart!')

        self.model_part_destination.ProcessInfo = model_part_origin.ProcessInfo # for detecting output writing

        for var in GenerateVariableListFromInput(settings["historical_variables_destination"]):
            self.model_part_destination.AddNodalSolutionStepVariable(var)

        self.data_comm = model_part_origin.GetCommunicator().GetDataCommunicator()

        self.destination_rank = settings["destination_rank"].GetInt()
        if self.destination_rank >= self.data_comm.Size():
            raise Exception("Destination rank %i larger than available size %i" %(self.destination_rank, self.data_comm.Size()))

        # optionally read mdpa (only on one rank)
        if mdpa_file_name_destination != "":
            if self.data_comm.Rank() == self.destination_rank:
                import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER
                KM.ModelPartIO(mdpa_file_name_destination, import_flags).ReadModelPart(self.model_part_destination)

            # properly initialize in MPI (on other ranks)
            if model_part_origin.IsDistributed():
                import KratosMultiphysics.mpi as KratosMPI

                # initialize SubModelPartStructure on other ranks
                KratosMPI.DistributedModelPartInitializer(self.model_part_destination, self.data_comm, self.destination_rank).CopySubModelPartStructure()

                data_comm_destination = KratosMPI.DataCommunicatorFactory.CreateFromRanksAndRegister(
                    self.data_comm,
                    [self.destination_rank],
                    "destination_mapping")

                if self.data_comm.Rank() != self.destination_rank:
                    KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicatorRecursively(self.model_part_destination, data_comm_destination)

        # optionally create output process (only on one rank)
        self.output_process = None
        if len(settings["output_process_settings"].keys()) > 0 and self.data_comm.Rank() == self.destination_rank:
            output_proc_params = KM.Parameters('''{ "dummy" : [] }''')
            output_proc_params["dummy"].Append(settings["output_process_settings"])
            self.output_process = KratosProcessFactory(model).ConstructListOfProcesses(output_proc_params["dummy"])[0]

        # create mapper
        if model_part_origin.IsDistributed():
            from KratosMultiphysics.MappingApplication.MPIExtension import MPIMapperFactory

            self.mapper = MPIMapperFactory.CreateMapper(
                model_part_origin,
                self.model_part_destination,
                settings["mapper_settings"])
        else:
            self.mapper = KM.MapperFactory.CreateMapper(
                model_part_origin,
                self.model_part_destination,
                settings["mapper_settings"])

    def ExecuteFinalizeSolutionStep(self):
        defaults = KM.Parameters('''{
            "variable_origin"      : "UNSPECIFIED",
            "variable_destination" : "UNSPECIFIED",
            "mapping_options" : []
        }''')

        for mapping_params in self.settings["mapping_settings"]:
            mapping_params.ValidateAndAssignDefaults(defaults)

            variable_origin = KM.KratosGlobals.GetVariable(mapping_params["variable_origin"].GetString())
            variable_destination = KM.KratosGlobals.GetVariable(mapping_params["variable_destination"].GetString())
            mapper_flags = GetMapperFlags(mapping_params["mapping_options"])

            self.mapper.Map(variable_origin, variable_destination, mapper_flags)

    def IsOutputStep(self):
        is_output_step = False
        if self.data_comm.Rank() == self.destination_rank:
            if self.output_process:
                is_output_step = self.output_process.IsOutputStep()
        return bool(self.data_comm.Broadcast(int(is_output_step), self.destination_rank))

    def PrintOutput(self):
        if self.data_comm.Rank() == self.destination_rank and self.output_process:
            self.output_process.PrintOutput()


def GetMapperFlags(settings):
    mapper_flags_dict = {
        "add_values"    : KM.Mapper.ADD_VALUES,
        "swap_sign"     : KM.Mapper.SWAP_SIGN,
        "use_transpose" : KM.Mapper.USE_TRANSPOSE
    }
    mapper_flags = KM.Flags()
    for flag_name in settings.GetStringArray():
        mapper_flags |= mapper_flags_dict[flag_name]

    return mapper_flags
