from pathlib import Path
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
import KratosMultiphysics.kratos_utilities as kratos_utils

def Factory(models, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    if not isinstance(parameters, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not ((isinstance(models, dict) and all(isinstance(x, Kratos.Model) for x in models.values())) or isinstance(models, Kratos.Model)):
        raise Exception("expected input shall be a dictionary of model objects or at least a model object alone")
    if not parameters.Has("model_name"):
        raise RuntimeError(f"VtuOutputProcess instantiation requires a \"model_name\" in parameters [ parameters = {parameters}].")
    if isinstance(models, dict):
        return VtuOutputProcess(models[parameters["model_name"].GetString()], parameters["settings"])
    else:
        return VtuOutputProcess(models, parameters["settings"])


class VtuOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""
        {
            "model_part_name"                   : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "file_format"                       : "binary",
            "output_precision"                  : 7,
            "output_control_type"               : "step",
            "output_interval"                   : 1.0,
            "output_sub_model_parts"            : false,
            "output_path"                       : {},
            "save_output_files_in_folder"       : true,
            "delete_existing_data_folder"       : false,
            "write_deformed_configuration"      : false,
            "nodal_solution_step_data_variables": [],
            "nodal_data_value_variables"        : [],
            "nodal_flags"                       : [],
            "element_data_value_variables"      : [],
            "element_flags"                     : [],
            "condition_data_value_variables"    : [],
            "condition_flags"                   : []

        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        super().__init__()

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.model_part_name = parameters["model_part_name"].GetString()
        self.write_deformed_configuration = parameters["write_deformed_configuration"].GetBool()
        self.output_precision = parameters["output_precision"].GetInt()
        self.output_sub_model_parts = parameters["output_sub_model_parts"].GetBool()
        self.nodal_solution_step_data_variables = parameters["nodal_solution_step_data_variables"].GetStringArray()
        self.nodal_data_value_variables = parameters["nodal_data_value_variables"].GetStringArray()
        self.nodal_flags = parameters["nodal_flags"].GetStringArray()
        self.element_data_value_variables = parameters["element_data_value_variables"].GetStringArray()
        self.element_flags = parameters["element_flags"].GetStringArray()
        self.condition_data_value_variables = parameters["condition_data_value_variables"].GetStringArray()
        self.condition_flags = parameters["condition_flags"].GetStringArray()

        file_format = parameters["file_format"].GetString()
        if file_format == "ascii":
            self.writer_format = Kratos.VtuOutput.ASCII
        elif file_format == "binary":
            self.writer_format = Kratos.VtuOutput.BINARY
        else:
            raise RuntimeError(f"Unsupported file format requested [ requested format = {file_format} ]. Supported file formats:\n\tascii\n\tbinary")

        if parameters["save_output_files_in_folder"].GetBool():
            self.output_path_training_samples = Path(parameters["output_path"]["train"].GetString())
            self.output_path_validation_samples = Path(parameters["output_path"]["valid"].GetString())
            if parameters["delete_existing_data_folder"].GetBool():
                if str(self.output_path_training_samples) != "":
                    kratos_utils.DeleteDirectoryIfExisting(str(self.output_path_training_samples))
                    if str(self.output_path_training_samples) != "":
                        self.output_path_training_samples.mkdir(parents=True, exist_ok=True)
                if str(self.output_path_validation_samples) != "":
                    kratos_utils.DeleteDirectoryIfExisting(str(self.output_path_validation_samples))
                    if str(self.output_path_validation_samples) != "":
                        self.output_path_validation_samples.mkdir(parents=True, exist_ok=True)
            #self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
            # now create the output path
            
        else:
            self.output_path = Path(".")

        self.vtu_output_ios: 'list[Kratos.VtuOutput]' = []

        #self.__controller = Kratos.OutputController(model, parameters)

    def ExecuteInitialize(self) -> None:
        # check and create all the vtu outputs
        self.model_part = self.model[self.model_part_name]
        self.vtu_output_ios.append(Kratos.VtuOutput(self.model_part, not self.write_deformed_configuration, self.writer_format, self.output_precision))

        if self.output_sub_model_parts:
            for sub_model_part in self.model_part.SubModelParts:
                self.vtu_output_ios.append(Kratos.VtuOutput(sub_model_part, not self.write_deformed_configuration, self.writer_format, self.output_precision))

        for vtu_output_io in self.vtu_output_ios:
            self.__AddData(vtu_output_io)
    
    def PrintOutput(self,type,sampleNo):
        '''type : training or validation
           sampleNo : nth sample '''
        for vtu_output in self.vtu_output_ios:
            vtu_output.PrintOutput(str(self.output_path / vtu_output.GetModelPart().FullName()) + "_" + type + "_{}".format(sampleNo))

    def PrintOutputTrainingSamples(self,sampleNo) -> None:
        ''' sampleNo : nth sample '''

        #current_suffix = self.__controller.GetCurrentControlValue()
        for vtu_output in self.vtu_output_ios:
            vtu_output.PrintOutput(str(self.output_path_training_samples / vtu_output.GetModelPart().FullName()) + "_{}".format(sampleNo))
        #self.__controller.Update()
    
    def PrintOutputValidationSamples(self,sampleNo) -> None:
        ''' sampleNo : nth sample '''

        #current_suffix = self.__controller.GetCurrentControlValue()
        for vtu_output in self.vtu_output_ios:
            vtu_output.PrintOutput(str(self.output_path_validation_samples / vtu_output.GetModelPart().FullName()) + "_{}".format(sampleNo))
        #self.__controller.Update()

    def IsOutputStep(self) -> bool:
        return self.__controller.Evaluate()

    def __AddData(self, vtu_output_io: Kratos.VtuOutput) -> None:
        for variable in self.nodal_solution_step_data_variables:
            vtu_output_io.AddHistoricalVariable(Kratos.KratosGlobals.GetVariable(variable))
        for variable in self.nodal_data_value_variables:
            vtu_output_io.AddNonHistoricalVariable(Kratos.KratosGlobals.GetVariable(variable), Kratos.VtuOutput.NODES)
        for variable in self.condition_data_value_variables:
            vtu_output_io.AddNonHistoricalVariable(Kratos.KratosGlobals.GetVariable(variable), Kratos.VtuOutput.CONDITIONS)
        for variable in self.element_data_value_variables:
            vtu_output_io.AddNonHistoricalVariable(Kratos.KratosGlobals.GetVariable(variable), Kratos.VtuOutput.ELEMENTS)

        for flag in self.nodal_flags:
            vtu_output_io.AddFlagVariable(flag, Kratos.KratosGlobals.GetFlag(flag), Kratos.VtuOutput.NODES)
        for flag in self.condition_flags:
            vtu_output_io.AddFlagVariable(flag, Kratos.KratosGlobals.GetFlag(flag), Kratos.VtuOutput.CONDITIONS)
        for flag in self.element_flags:
            vtu_output_io.AddFlagVariable(flag, Kratos.KratosGlobals.GetFlag(flag), Kratos.VtuOutput.ELEMENTS)

