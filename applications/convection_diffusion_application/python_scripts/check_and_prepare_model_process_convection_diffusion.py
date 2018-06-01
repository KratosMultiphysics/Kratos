from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    """Prepare the computing model part.

    The computing model part is created if it does not exist. Nodes and elements
    from the domain sub model parts are added to the computing model part.
    Conditions are added from the processes sub model parts.
    """
    def __init__(self, main_model_part, Parameters):
        self.main_model_part = main_model_part
        
        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        self.problem_domain_sub_model_part_list = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

    def Execute(self):
        
        problem_domain_sub_model_parts = []
        for i in range(self.problem_domain_sub_model_part_list.size()):
            problem_domain_sub_model_parts.append(self.main_model_part.GetSubModelPart(self.problem_domain_sub_model_part_list[i].GetString()))
        
        processes_parts = []
        for i in range(self.processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        #construct a model part which contains both the skin and the volume
        if (self.main_model_part.HasSubModelPart(self.computing_model_part_name)):
            convection_diffusion_computational_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        else:
            convection_diffusion_computational_model_part = self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        convection_diffusion_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        convection_diffusion_computational_model_part.Properties  = self.main_model_part.Properties
        
        for part in problem_domain_sub_model_parts:
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(convection_diffusion_computational_model_part, part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
            transfer_process.Execute()
        for part in processes_parts:
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(convection_diffusion_computational_model_part, part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS)
            transfer_process.Execute()

        KratosMultiphysics.Logger.PrintInfo("Computing model part:", convection_diffusion_computational_model_part)
