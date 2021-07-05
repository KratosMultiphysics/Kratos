import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.masking_process import MaskingProcess


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MaskZerosProcess(settings["parameters"])

class MaskZerosProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """

        self.objective = settings["objective"].GetString()
        try:
            self.log_denominator = settings["log_denominator"].GetString()
        except RuntimeError:
            self.log_denominator = "masking_zeros"
        
    def Preprocess(self, data_in, data_out):

        mask_parameters = KM.Parameters()
        # Setup for input
        if self.objective == 'input':
            dimension_variables = range(len(data_in[0]))
            dataset_length = range(len(data_in[:,0]))
            masking_variables = []
            # Extract which variables are always 0
            for variable in dimension_variables:
                if all(abs(data_in[i,variable]) < 1e-13 for i in dataset_length):
                    masking_variables.append(variable)
            try:
                mask_parameters.AddEmptyValue("input_log_name")
                mask_parameters["input_log_name"].SetString(self.input_log_name)
            except AttributeError:
                print("No input_log_name defined.")
                mask_parameters.RemoveValue("input_log_name")

        # Setup for output
        elif self.objective == 'output':
            dimension_variables = range(len(data_out[0]))
            dataset_length = range(len(data_out[:,0]))
            masking_variables = []
            # Extract which variables are always 0
            for variable in dimension_variables:
                if all(abs(data_out[i,variable]) < 1e-13 for i in dataset_length):
                    masking_variables.append(variable)
            try:
                mask_parameters.AddEmptyValue("output_log_name")
                mask_parameters["output_log_name"].SetString(self.output_log_name)
            except AttributeError:
                    print("No output_log_name defined.")
                    mask_parameters.RemoveValue("output_log_name")
        else:
            raise Exception("Masking objective not supported. Supported objectives are input and output")

        if len(masking_variables) > 0:
            # Load the parameters to the masking process
            mask_parameters.AddEmptyValue("objective")
            mask_parameters["objective"].SetString(self.objective)
            mask_parameters.AddEmptyValue("load_from_log")
            mask_parameters["load_from_log"].SetBool(self.load_from_log)
            mask_parameters.AddEmptyValue("variable_ids")
            mask_parameters["variable_ids"].SetVector(masking_variables)
            mask_parameters.AddEmptyValue("log_denominator")
            mask_parameters["log_denominator"].SetString(self.log_denominator)

            # Mask the dataset
            self.mask_process = MaskingProcess(mask_parameters)
            [data_in , data_out] = self.mask_process.Preprocess(data_in, data_out)

        return [data_in, data_out]
