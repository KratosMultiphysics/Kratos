from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics 

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    def __init__(self, main_model_part, Parameters ):
        self.main_model_part = main_model_part
        
        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        self.sub_model_part_names       = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

        self.bodies_parts_list = []
        self.bodies_list = False
        if( Parameters.Has("bodies_list") ):
            self.bodies_list = True
            self.bodies_parts_list = Parameters["bodies_list"]
        
        
    def Execute(self):

        #construct body model parts:
        solid_body_model_parts = []
        fluid_body_model_parts = []
        rigid_body_model_parts = []

        #construct body model parts:
        if( self.bodies_list == True ):
            for i in range(self.bodies_parts_list.size()):
                #create body model part
                body_model_part_name = self.bodies_parts_list[i]["body_name"].GetString()
                self.main_model_part.CreateSubModelPart(body_model_part_name)
                body_model_part = self.main_model_part.GetSubModelPart(body_model_part_name)
                
                print("::[Model_Prepare]::Body Created :", body_model_part_name)
                body_model_part.ProcessInfo = self.main_model_part.ProcessInfo
                body_model_part.Properties  = self.main_model_part.Properties

                #build body from their parts
                body_parts_name_list = self.bodies_parts_list[i]["parts_list"]
                body_parts_list = []
                #for j in range(len(body_parts_name_list)):
                for j in range(body_parts_name_list.size()):

                    body_parts_list.append(self.main_model_part.GetSubModelPart(body_parts_name_list[j].GetString()))

                body_model_part_type = self.bodies_parts_list[i]["body_type"].GetString()

                for part in body_parts_list:
                    for node in part.Nodes:
                        if (body_model_part_type=="Solid"):
                            node.Set(KratosMultiphysics.SOLID)
                        if (body_model_part_type=="Fluid"):
                            node.Set(KratosMultiphysics.FLUID)
                        if (body_model_part_type=="Rigid"):
                            node.Set(KratosMultiphysics.RIGID)
                
                for part in body_parts_list:
                    for node in part.Nodes:
                        body_model_part.Nodes.append(node)
                    for elem in part.Elements:
                        body_model_part.AddElement(elem,0)
                    for cond in part.Conditions:
                        body_model_part.AddCondition(cond,0) 

                if( body_model_part_type == "Solid" ):
                    body_model_part.Set(KratosMultiphysics.SOLID)
                    solid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
                if( body_model_part_type == "Fluid" ):
                    body_model_part.Set(KratosMultiphysics.FLUID)
                    fluid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
                if( body_model_part_type == "Rigid" ):
                    body_model_part.Set(KratosMultiphysics.RIGID)
                    rigid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))

            #add walls in fluid domains:
            for fluid_part in fluid_body_model_parts:
                for rigid_part in rigid_body_model_parts:
                    for node in rigid_part.Nodes:
                        if( node.IsNot(KratosMultiphysics.FLUID) ):
                            fluid_part.AddNode(node,0)
                            
        #construct the computing model part:
        domain_parts = []
        for i in range(self.sub_model_part_names.size()):
            domain_parts.append(self.main_model_part.GetSubModelPart(self.sub_model_part_names[i].GetString()))
        processes_parts = []
        for i in range(self.processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        #construct a model part which contains the mesh to compute
        self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        computing_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        computing_model_part.Properties  = self.main_model_part.Properties

        #set flag to identify the solid model part :: solid application
        computing_model_part.Set(KratosMultiphysics.SOLID)
        #set flag to identify the computing model part
        computing_model_part.Set(KratosMultiphysics.ACTIVE)
        
        for node in self.main_model_part.Nodes:
            computing_model_part.AddNode(node,0)
 
        for part in domain_parts:
            #part.Set(KratosMultiphysics.SOLID)
            for elem in part.Elements:
                computing_model_part.AddElement(elem,0)
            #for node in part.Nodes:
            #    computing_model_part.Nodes.append(node,0)
            
        for part in processes_parts:
            part.Set(KratosMultiphysics.BOUNDARY)
            for cond in part.Conditions:
                computing_model_part.AddCondition(cond,0)


        #delete body parts: (materials have to be already assigned)
        if( self.bodies_list == True ):
            for i in range(self.bodies_parts_list.size()):
               #get body parts
                body_parts_name_list = self.bodies_parts_list[i]["parts_list"]
                for j in range(body_parts_name_list.size()):
                    self.main_model_part.RemoveSubModelPart(body_parts_name_list[j].GetString())
                    print("::[Model_Prepare]::Body Part Removed:", body_parts_name_list[j].GetString())
       
        #for part in domain_parts:
        #    self.main_model_part.RemoveSubModelPart(part)
        #    print("Removed SubModelPart:", part.Name)
          
        print("::[Model_Prepare]::",computing_model_part)       
