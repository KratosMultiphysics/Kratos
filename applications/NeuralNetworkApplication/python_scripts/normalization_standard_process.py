import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.centering_process import CenteringProcess
from KratosMultiphysics.NeuralNetworkApplication.scaling_process import ScalingProcess


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return NormalizationStandardProcess(settings["parameters"])

class NormalizationStandardProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing process settings.
        """

        self.objective = settings["objective"].GetString()
        try:
            self.log_denominator = settings["log_denominator"].GetString()
        except RuntimeError:
            self.log_denominator = "normalization"
        
        # Centering
        if settings.Has("center"):
            self.center = settings["center"].GetString()
        else:
            self.center = "mean"
        centering_parameters = KM.Parameters()
        centering_parameters.AddEmptyValue("center")
        centering_parameters["center"].SetString(self.center)
        centering_parameters.AddEmptyValue("objective")
        centering_parameters["objective"].SetString(self.objective)
        try:
            centering_parameters.AddEmptyValue("input_log_name")
            centering_parameters["input_log_name"].SetString(self.input_log_name)
            centering_parameters.AddEmptyValue("output_log_name")
            centering_parameters["output_log_name"].SetString(self.output_log_name)
        except AttributeError:
            centering_parameters.RemoveValue("input_log_name")
            centering_parameters.RemoveValue("output_log_name")
        centering_parameters.AddEmptyValue("load_from_log")
        centering_parameters["load_from_log"].SetBool(self.load_from_log)
        centering_parameters.AddEmptyValue("log_denominator")
        centering_parameters["log_denominator"].SetString(self.log_denominator + 'center')
        self.center_process = CenteringProcess(centering_parameters)

        # Scaling
        if settings.Has("scale"):
            self.scale = settings["scale"].GetString()
        else:
            self.scale = "std"
        scaling_parameters = KM.Parameters()
        scaling_parameters.AddEmptyValue("scale")
        scaling_parameters["scale"].SetString(self.scale)
        scaling_parameters.AddEmptyValue("objective")
        scaling_parameters["objective"].SetString(self.objective)
        try:
            scaling_parameters.AddEmptyValue("input_log_name")
            scaling_parameters["input_log_name"].SetString(self.input_log_name)
            scaling_parameters.AddEmptyValue("output_log_name")
            scaling_parameters["output_log_name"].SetString(self.output_log_name)
        except AttributeError:
            scaling_parameters.RemoveValue("input_log_name")
            scaling_parameters.RemoveValue("output_log_name")
        scaling_parameters.AddEmptyValue("load_from_log")
        scaling_parameters["load_from_log"].SetBool(self.load_from_log)
        scaling_parameters.AddEmptyValue("log_denominator")
        scaling_parameters["log_denominator"].SetString(self.log_denominator + 'scaling')
        self.scale_process = ScalingProcess(scaling_parameters)
        
    def Preprocess(self, data_in, data_out):
        
        [data_in , data_out] = self.center_process.Preprocess(data_in, data_out)
        [data_in , data_out] = self.scale_process.Preprocess(data_in, data_out)

        return [data_in, data_out]
    
    def Invert(self, data_in, data_out):
        
        [data_in , data_out] = self.scale_process.Invert(data_in, data_out)
        [data_in , data_out] = self.center_process.Invert(data_in, data_out)

        return [data_in, data_out]

