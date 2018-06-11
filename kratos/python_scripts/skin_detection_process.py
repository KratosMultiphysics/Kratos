from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics 
 
def Factory(settings, Model): 
    if(type(settings) != KratosMultiphysics.Parameters): 
        raise Exception("expected input shall be a Parameters object, encapsulating a json string") 
    return SkinDetectionProcess(Model, settings["Parameters"])
 
##all the processes python processes should be derived from "python_process" 
class SkinDetectionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ): 
        KratosMultiphysics.Process.__init__(self) 
 
        default_settings = KratosMultiphysics.Parameters(""" 
        {
            "model_part_name"          : "Main",
            "recursive_detection"      : false,
            "name_auxiliar_model_part" : "SkinModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 0
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.recursive_detection = settings["recursive_detection"].GetBool()

        detect_skin_parameters = KratosMultiphysics.Parameters("""{}""")
        detect_skin_parameters.AddValue("name_auxiliar_model_part", settings["name_auxiliar_model_part"])
        detect_skin_parameters.AddValue("name_auxiliar_condition", settings["name_auxiliar_condition"])
        detect_skin_parameters.AddValue("echo_level", settings["echo_level"])
        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess2D(self.model_part, detect_skin_parameters)
        else:
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess3D(self.model_part, detect_skin_parameters)

    def ExecuteInitialize(self):
        # We execute the process
        self.detect_skin.Execute()

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        # If recurive we detect each time step
        if (self.recursive_detection is True):
            self.detect_skin.Execute()

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def Clear(self):
        pass
