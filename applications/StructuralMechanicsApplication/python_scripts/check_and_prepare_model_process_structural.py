from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

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
        self.structural_model_part_names = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

    def Execute(self):
        
        structural_parts = []
        for i in range(self.structural_model_part_names.size()):
            structural_parts.append(self.main_model_part.GetSubModelPart(self.structural_model_part_names[i].GetString()))
        
        processes_parts = []
        for i in range(self.processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        #construct a model part which contains both the skin and the volume
        if (self.main_model_part.HasSubModelPart(self.computing_model_part_name)):
            structural_computational_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        else:
            structural_computational_model_part = self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        structural_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        structural_computational_model_part.Properties  = self.main_model_part.Properties
        
        #node_ids = set()
        #elem_ids = set()
        #for part in structural_parts:
            #for node in part.Nodes:
                #node_ids.add(node.Id)
            #for elem in part.Elements:
                #elem_ids.add(elem.Id)
                
        #structural_computational_model_part.AddNodes(list(node_ids))
        #structural_computational_model_part.AddElements(list(elem_ids))
            
        #cond_ids = set()
        #for part in processes_parts:
            #for cond in part.Conditions:
                #cond_ids.add(cond.Id)
        #structural_computational_model_part.AddConditions(list(cond_ids)) 
        
        for part in structural_parts:
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(structural_computational_model_part, part, "NodesAndElements")
            transfer_process.Execute()
        for part in processes_parts:
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(structural_computational_model_part, part, "Conditions")
            transfer_process.Execute()

        print("Computing model part:")        
        print(structural_computational_model_part)
