# Importing the Kratos Library
import KratosMultiphysics as KM

def Factory(settings, Model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class CheckAndPrepareModelProcess(KM.Process):
    """Prepare the computing model part.

    The computing model part is created if it does not exist. Nodes and elements
    from the domain sub model parts are added to the computing model part.
    Conditions are added from the processes sub model parts.
    """
    def __init__(self, Model, Parameters):
        KM.Process.__init__(self)
        self.main_model_part_name = Parameters["model_part_name"].GetString()
        self.main_model_part = Model[self.main_model_part_name]
        self.model = Model

        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        self.structural_model_part_names = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

    def Execute(self):

        structural_parts = self.__ExtractModelparts(self.structural_model_part_names)
        processes_parts  = self.__ExtractModelparts(self.processes_model_part_names)

        # Construct a model part which contains both the skin and the volume
        if (self.main_model_part.HasSubModelPart(self.computing_model_part_name)):
            structural_computational_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        else:
            structural_computational_model_part = self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        structural_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        structural_computational_model_part.Properties  = self.main_model_part.Properties

        for part in structural_parts:
            KM.FastTransferBetweenModelPartsProcess(structural_computational_model_part, part,
                KM.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS).Execute()

        for part in processes_parts:
            KM.FastTransferBetweenModelPartsProcess(structural_computational_model_part, part,
                KM.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS).Execute()

        KM.Logger.PrintInfo("Computing model part", structural_computational_model_part)

    def __ExtractModelparts(self, params_list):
        """Extracting a list of SubModelParts to be added to the ComputingModelPart
        If the MainModelPart is to be added to the ComputingModelPart, then this is
        done directly, no need to also add the SubModelParts, since they are contained
        in the MainModelpart
        """
        list_model_part_names = [params_list[i].GetString() for i in range(params_list.size())]
        if self.main_model_part_name in list_model_part_names:
            # directly return if the MainModelPart is to be added to the ComputingModelPart
            return [self.main_model_part]
        else:
            model_parts_list = []
            for model_part_name  in list_model_part_names:
                # only SubModelParts of the MainModelPart can be used!
                full_model_part_name = self.main_model_part_name + "." + model_part_name
                model_parts_list.append(self.model[full_model_part_name])
            return model_parts_list
