# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return FlowOutputProcess(Model, settings["Parameters"])

class FlowOutputProcess(KratosMultiphysics.Process):

    """This process calculates(using c++ utilities) and writes the flow through a given list of (sub)model parts.
    In 3D use a surface eg. Inlet, Outlet
    In 2D use a line eg. Inlet, Outlet

    This process works in MPI as well as with restarts
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters('''{
            "help"            : "This process calculates(using c++ utilities) and writes the flow through a given list of (sub)model parts.",
            "model_part_name_list" : [],
            "print_format"      : "",
            "output_file_settings": {}
        }''')

        self.model = model
        self.params = params
        self.params.ValidateAndAssignDefaults(default_settings)
        self.output_file = None
        self.format = self.params["print_format"].GetString()

    def ExecuteInitialize(self):
        # getting the ModelPart from the Model
        model_part_name_list = self.params["model_part_name_list"]
        if model_part_name_list.size() == 0:
            raise Exception('No model parts are specified!')

        self.model_part_for_time = self.model[model_part_name_list[0].GetString()]

        # Only rank 0 writes in MPI
        my_rank = 0
        comm = self.model_part_for_time.GetCommunicator().GetDataCommunicator()
        self.is_writing_rank = my_rank == comm.Rank()
        if self.is_writing_rank:
            file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])
            file_header = self.GetFileHeader()
            self.output_file =  TimeBasedAsciiFileWriterUtility(self.model_part_for_time, file_handler_params, file_header).file

    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part_for_time.ProcessInfo[KratosMultiphysics.TIME]
        model_part_name_list = self.params["model_part_name_list"]

        out = str(time)
        for model_part_name_param in model_part_name_list:
            model_part_name = model_part_name_param.GetString()
            model_part = self.model[model_part_name]
            flow_value = self.CalculateFlow(model_part)

            out += " " + format(flow_value,self.format)

        out += "\n"
        if self.is_writing_rank:
            self.output_file.write(out)

    def ExecuteFinalize(self):
        if self.is_writing_rank:
            self.output_file.close()

    def GetFileHeader(self):
        model_part_name_list = self.params["model_part_name_list"]
        header  = '# Flow results ' + '\n'
        header += '# time '
        for model_part_name_param in model_part_name_list:
            model_part_name = model_part_name_param.GetString()
            model_part = self.model[model_part_name]
            header += model_part.Name
            header += ' '
        header += "\n"

        return header

    def CalculateFlow(self, model_part):
        flow_value = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRate(model_part)
        return flow_value
