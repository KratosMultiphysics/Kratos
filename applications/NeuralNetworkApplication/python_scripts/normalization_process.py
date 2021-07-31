import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.normalization_standard_process import NormalizationStandardProcess
from KratosMultiphysics.NeuralNetworkApplication.normalization_zero_one_process import NormalizationZeroOneProcess
from KratosMultiphysics.NeuralNetworkApplication.normalization_minus_one_one_process import NormalizationMinusOneOneProcess
from KratosMultiphysics.NeuralNetworkApplication.normalization_soft_zero_one_process import NormalizationSoftZeroOneProcess


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return NormalizationProcess(settings["parameters"])

class NormalizationProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing process settings.
        """

        self.objective = settings["objective"].GetString()
        try:
            self.normalization_type = settings["normalization_type"].GetString()
        except RuntimeError:
            self.normalization_type = "standard"
        try:
            self.log_denominator = settings["log_denominator"].GetString()
        except RuntimeError:
            self.log_denominator = "normalization"
        
        normalization_parameters = KM.Parameters()
        normalization_parameters.AddEmptyValue("objective")
        normalization_parameters["objective"].SetString(self.objective)
        try:
            normalization_parameters.AddEmptyValue("input_log_name")
            normalization_parameters["input_log_name"].SetString(self.input_log_name)
            normalization_parameters.AddEmptyValue("output_log_name")
            normalization_parameters["output_log_name"].SetString(self.output_log_name)
        except AttributeError:
            normalization_parameters.RemoveValue("input_log_name")
            normalization_parameters.RemoveValue("output_log_name")
        normalization_parameters.AddEmptyValue("load_from_log")
        normalization_parameters["load_from_log"].SetBool(self.load_from_log)
        normalization_parameters.AddEmptyValue("log_denominator")
        normalization_parameters["log_denominator"].SetString(self.log_denominator)

        # Centering
        if settings.Has("center"):
            self.center = settings["center"].GetString()
            normalization_parameters.AddEmptyValue("center")
            normalization_parameters["center"].SetString(self.center)
        # Centering
        if settings.Has("scale"):
            self.scale = settings["scale"].GetString()
            normalization_parameters.AddEmptyValue("scale")
            normalization_parameters["scale"].SetString(self.scale)
            
        # Select normalization type
        switcher = {
            'standard': NormalizationStandardProcess,
            'zero_one': NormalizationZeroOneProcess,
            'minus_one_one': NormalizationMinusOneOneProcess,
            'soft_zero_one': NormalizationSoftZeroOneProcess
        }
        
        def __GetNormalization(norm_type, parameters):
            norm = switcher.get(norm_type, NormalizationStandardProcess)
            return norm(parameters)

        self.normalization_process = __GetNormalization(self.normalization_type, normalization_parameters)
        
    def Preprocess(self, data_in, data_out):
        
        [data_in , data_out] = self.normalization_process.Preprocess(data_in, data_out)

        return [data_in, data_out]
    
    def Invert(self, data_in, data_out):
        
        [data_in , data_out] = self.normalization_process.Invert(data_in, data_out)

        return [data_in, data_out]

