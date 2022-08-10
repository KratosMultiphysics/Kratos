# Importing the Kratos Library
from urllib import response
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, Model):
    if(type(settings) != Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ResponseFunctionOutputProcess(Model, settings["Parameters"])

class ResponseFunctionOutputProcess(Kratos.OutputProcess):
    def __init__(self, model, params):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters('''{
            "response_type"       : "PLEASE_SPECIFY_RESPONSE_TYPE",
            "model_part_name"     : "PLEASE_SPECIFY_MAIN_MODEL_PART_NAME",
            "response_settings"   : {},
            "output_file_settings": {}
        }''')

        self.model = model
        self.params = params
        self.main_model_part = self.model.GetModelPart(
            self.params["model_part_name"].GetString())
        self.params.ValidateAndAssignDefaults(default_settings)
        self.output_file = None
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

        response_type = self.params["response_type"].GetString()
        if (response_type == "norm_square"):
            self.response = KratosCFD.VelocityPressureNormSquareResponseFunction(
                self.params["response_settings"], self.model)
        elif (response_type == "domain_integrated"):
            self.response = KratosCFD.DomainIntegratedResponseFunction(
                self.params["response_settings"], self.main_model_part)
        elif (response_type == "domain_integrated_3d_vector_magnitude_square_power_mean"):
            if domain_size == 2:
                self.response = KratosCFD.DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction2D(
                    self.params["response_settings"], self.main_model_part)
            elif domain_size == 3:
                self.response = KratosCFD.DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction3D(
                    self.params["response_settings"], self.main_model_part)
            else:
                raise Exception("Unsupported domain size requested.")
        else:
            raise Exception(
                "Unknown response_type = \"" + response_type +
                "\". Supported response types are: \n\t  1. norm_square")

    def ExecuteInitialize(self):
        self.response.Initialize()
        # Only rank 0 writes in MPI
        my_rank = 0
        comm = self.main_model_part.GetCommunicator().GetDataCommunicator()
        self.is_writing_rank = my_rank == comm.Rank()
        if self.is_writing_rank:
            file_handler_params = Kratos.Parameters(self.params["output_file_settings"])
            file_header = self.GetFileHeader()
            self.output_file =  TimeBasedAsciiFileWriterUtility(self.main_model_part, file_handler_params, file_header).file

    def PrintOutput(self):
        time = self.main_model_part.ProcessInfo[Kratos.TIME]

        out = str(time)
        response_value = self.response.CalculateValue(self.main_model_part)
        out += "," + str(response_value)
        out += "\n"

        if self.is_writing_rank:
            self.output_file.write(out)

    def ExecuteFinalize(self):
        if self.is_writing_rank:
            self.output_file.close()

    def GetFileHeader(self):
        header  = '# Response results ' + '\n'
        header += '# time, ResponseValue'
        header += "\n"

        return header