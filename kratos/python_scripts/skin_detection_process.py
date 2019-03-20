from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics 
 
def Factory(settings, Model): 
    if(type(settings) != KratosMultiphysics.Parameters): 
        raise Exception("expected input shall be a Parameters object, encapsulating a json string") 
    return SkinDetectionProcess(Model, settings["Parameters"])
 
## All the processes python should be derived from "Process" 
class SkinDetectionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ): 
        KratosMultiphysics.Process.__init__(self) 
 
        default_settings = KratosMultiphysics.Parameters(""" 
        {
            "help"                                  : "This process detects the skin from a given submodelpart and it generates the correspong conditions",
            "model_part_name"                       : "Main",
            "computing_model_part_name"             : "computing_domain",
            "recursive_detection"                   : false,
            "name_auxiliar_model_part"              : "SkinModelPart",
            "name_auxiliar_condition"               : "Condition",
            "list_model_parts_to_assign_conditions" : [],
            "echo_level"                            : 0
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        # The main model part
        self.model_part = Model[settings["model_part_name"].GetString()]
        # The computing model part
        self.computing_model_part = Model[settings["computing_model_part_name"].GetString()]

        # Recursive detection
        self.recursive_detection = settings["recursive_detection"].GetBool()

        # Skin model part name
        self.name_auxiliar_model_part = settings["name_auxiliar_model_part"].GetString()

        # Process parameters
        detect_skin_parameters = KratosMultiphysics.Parameters("""{}""")
        detect_skin_parameters.AddValue("name_auxiliar_model_part", settings["name_auxiliar_model_part"])
        detect_skin_parameters.AddValue("name_auxiliar_condition", settings["name_auxiliar_condition"])
        detect_skin_parameters.AddValue("list_model_parts_to_assign_conditions", settings["list_model_parts_to_assign_conditions"])
        detect_skin_parameters.AddValue("echo_level", settings["echo_level"])
        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess2D(self.model_part, detect_skin_parameters)
        else:
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess3D(self.model_part, detect_skin_parameters)

    def ExecuteInitialize(self):
        # We execute the process
        self.detect_skin.Execute()
        # We copy the conditions to the contact model part
        skin_model_part = self.model_part.GetSubModelPart(self.name_auxiliar_model_part)
        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.computing_model_part, skin_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS)
        transfer_process.Execute()

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        # If recurive we detect each time step
        if (self.recursive_detection is True):
            self.detect_skin.Execute()
            # We copy the conditions to the contact model part
            skin_model_part = self.model_part.GetSubModelPart(self.name_auxiliar_model_part)
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.computing_model_part, skin_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS)
            transfer_process.Execute()

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
