# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMWriteEnergyOutputProcess(model, settings["Parameters"])

class MPMWriteEnergyOutputProcess(KratosMultiphysics.OutputProcess):
    """This process writes the (kinetic, potential, strain and total) energy of a model part.
    """

    def __init__(self, model, params):
        super().__init__()

        self.model = model
        self.params = params
        self.params.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.interval = KratosMultiphysics.IntervalUtility(params)
        self.output_file = None
        self.format = self.params["print_format"].GetString()

        # Getting the ModelPart from the Model
        self.model_part_name = self.params["model_part_name"].GetString()
        if not self.model_part_name:
            raise Exception('No "model_part_name" was specified!')
        elif self.model_part_name.startswith('Background_Grid'):
            self.model_part_name = self.full_model_part_name.replace('Background_Grid','MPM_Material')

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters('''{
            "model_part_name"      : "",
            "interval"             : [0.0, 1e30],
            "print_format"         : ".8f",
            "output_file_settings" : {}
        }''')

    def ExecuteBeforeSolutionLoop(self):
        self.model_part = self.model[self.model_part_name]
        file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])

        # if output file name is not specified, use default value: model_part_name_energy.dat
        if not file_handler_params.Has("file_name"):
            output_file_name = f"{self.model_part_name}_energy.dat"
            wrn_msg = f"File name not specified, using the default name: {output_file_name}"
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__,wrn_msg)
            file_handler_params.AddString("file_name",output_file_name)

        file_header = self._GetFileHeader()
        self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_params, file_header).file
        if self.output_file:
            self.output_file.flush()

    def IsOutputStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        return self.interval.IsInInterval(time)

    def PrintOutput(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if self.output_file:
            p_energy = k_energy = s_energy = t_energy = 0.0
            KratosMPM.EnergyCalculationUtility().CalculateTotalEnergy(self.model_part,p_energy,k_energy,s_energy,t_energy)
            out = f"{str(time)} {p_energy:{self.format}} {k_energy:{self.format}} {s_energy:{self.format}} {t_energy:{self.format}}\n"
            self.output_file.write(out)
            self.output_file.flush()

    def ExecuteFinalize(self):
        if self.output_file:
            self.output_file.close()

    def _GetFileHeader(self):
        header  = f"# Energy model part '{self.model_part_name}'\n"
        header += f"# time potential_energy kinetic_energy strain_energy total_energy\n"
        return header
