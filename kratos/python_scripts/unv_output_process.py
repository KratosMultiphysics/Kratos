import KratosMultiphysics
from KratosMultiphysics import kratos_utilities
from  KratosMultiphysics.deprecation_management import DeprecationManager

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model):
    if(isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if(isinstance(model, KratosMultiphysics.Model)):
        raise Exception("expected input shall be a Model object")

    return UnvOutputProcess(model, settings["Parameters"])

class UnvOutputProcess(KratosMultiphysics.Process):
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters ):
        KratosMultiphysics.Process.__init__(self)

        # IMPORTANT: when "output_control_type" is "time",
        # then paraview will not be able to group them
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_interval"                    : 1.0,
            "output_control_type"                : "step",
            "nodal_results"                      : []
        }""")

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        self.settings = settings

        # Warning: we may be changing the parameters object here:
        self.TranslateLegacyVariablesAccordingToCurrentStandard(settings)

        self.settings.ValidateAndAssignDefaults(default_parameters)

        # if self.settings["save_output_files_in_folder"].GetBool():
        #     if self.model_part.GetCommunicator().MyPID() == 0:
        #         folder_name = self.settings["folder_name"].GetString()
        #         if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
        #             import KratosMultiphysics.kratos_utilities as kratos_utils
        #             kratos_utils.DeleteDirectoryIfExisting(folder_name)
        #         if not os.path.isdir(folder_name):
        #             os.mkdir(folder_name)
        #     self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        # self.unv_io = KratosMultiphysics.UnvOutput(self.model_part, self.settings)
        self.unv_io = KratosMultiphysics.UnvOutput(self.model_part, "nxout")
        self.unv_io.InitializeMesh()
        self.unv_io.WriteMesh()

        output_control_type = settings["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_variable = KratosMultiphysics.TIME
            self.output_control_utility = KratosMultiphysics.DoubleFixedIntervalRecurringEventUtility(self.model_part.ProcessInfo[self.output_control_variable], settings["output_interval"].GetDouble())
        elif output_control_type == "step":
            self.output_control_variable = KratosMultiphysics.STEP
            self.output_control_utility = KratosMultiphysics.IntegerFixedIntervalRecurringEventUtility(self.model_part.ProcessInfo[self.output_control_variable], settings["output_interval"].GetInt())
        else:
            raise RuntimeError(f"Unsupported output control type = \"{output_control_type}\" requested. Supported control types are:\n\ttime\n\tstep")

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings):

        context_string = type(self).__name__

        old_name = 'output_frequency'
        new_name = 'output_interval'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

    def ExecuteInitialize(self):
        self.nodal_variables = kratos_utilities.GenerateVariableListFromInput(self.settings["nodal_results"])

    def PrintOutput(self):
        current_control_value = self.model_part.ProcessInfo[self.output_control_variable]
        for variable in self.nodal_variables:
            self.unv_io.PrintOutput(variable, current_control_value)
        self.output_control_utility.ScheduleNextEvent(current_control_value)

    def IsOutputStep(self):
        return self.output_control_utility.IsEventExpected(self.model_part.ProcessInfo[self.output_control_variable])