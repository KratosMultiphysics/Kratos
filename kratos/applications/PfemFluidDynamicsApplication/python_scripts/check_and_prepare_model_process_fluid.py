#from KratosMultiphysics import *
import KratosMultiphysics 
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemBaseApplication import *
from KratosMultiphysics.PfemFluidDynamicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
#from KratosMultiphysics.PfemBaseApplication import *
#from KratosMultiphysics.PfemFluidDynamicsApplication import *

import time as timer

def StartTimeMeasuring():
    # Measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # Measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[PFEM_FLUID_MODEL]:: [ %.2f" % round(used_time,2),"s", process," ] ")


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    def __init__(self, main_model_part, Parameters ):
        self.main_model_part = main_model_part
        
        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        if( Parameters.Has("problem_domain_sub_model_part_list") ):
            self.sub_model_part_names     = Parameters["problem_domain_sub_model_part_list"]
        if( Parameters.Has("processes_sub_model_part_list") ):
            self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

        self.bodies_parts_list = []
        self.bodies_list = False
        if( Parameters.Has("bodies_list") ):
            self.bodies_list = True
            self.bodies_parts_list = Parameters["bodies_list"]


    def Execute(self):
        print("Execute check_and_prepare_model_process_fluid")


        #construct body model parts:
        fluid_body_model_parts = []
        solid_body_model_parts = []
        rigid_body_model_parts = []

        if( self.bodies_list == True ):
            for i in range(self.bodies_parts_list.size()):
                #create body model part
                body_model_part_name = self.bodies_parts_list[i]["body_name"].GetString()
                self.main_model_part.CreateSubModelPart(body_model_part_name)
                body_model_part = self.main_model_part.GetSubModelPart(body_model_part_name)
 
                print("[Model_Prepare]::Body Creation", body_model_part_name)
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
                    
                    clock_time = StartTimeMeasuring()
                    for node in part.Nodes:
                        #body_model_part.AddNode(node,0)
                        body_model_part.Nodes.append(node)
                        if (body_model_part_type=="Fluid"):
                            node.Set(KratosMultiphysics.FLUID)
                        if (body_model_part_type=="Solid"):
                            node.Set(KratosMultiphysics.FLUID)
                            node.Set(KratosMultiphysics.SOLID)
                        if (body_model_part_type=="Rigid"):
                            node.Set(KratosMultiphysics.RIGID)
                            node.Set(KratosMultiphysics.BOUNDARY)
                    StopTimeMeasuring(clock_time,"1. for node in part.Nodes", True);

                    clock_time = StartTimeMeasuring()
                    for elem in part.Elements:
                        #body_model_part.AddElement(elem,0)
                        body_model_part.Elements.append(elem)
                    StopTimeMeasuring(clock_time,"1. part.Elements", True);
                    clock_time = StartTimeMeasuring()
                    for cond in part.Conditions:
                        #body_model_part.AddCondition(cond,0) 
                        body_model_part.Conditions.append(cond) 
                    StopTimeMeasuring(clock_time,"1. part.Conditions", True);

                if( body_model_part_type == "Fluid" ):
                    body_model_part.Set(KratosMultiphysics.FLUID)
                    fluid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
                if( body_model_part_type == "Solid" ):
                    body_model_part.Set(KratosMultiphysics.FLUID)
                    body_model_part.Set(KratosMultiphysics.SOLID)
                    solid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
                if( body_model_part_type == "Rigid" ):
                    body_model_part.Set(KratosMultiphysics.RIGID)
                    rigid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))

                print(body_model_part)

            clock_time = StartTimeMeasuring()

            #for fluid_part in fluid_body_model_parts:
            #    for rigid_part in rigid_body_model_parts:
            #        for node in rigid_part.Nodes:
            #            if( node.IsNot(KratosMultiphysics.FLUID) ):
            #                #fluid_part.AddNode(node,0)
            #                fluid_part.Nodes.append(node)
            #                print("Node Inserted Py",node.Id)

            #add walls in fluid domains:
            node_flags = FlagsContainer()
            #node_flags.PushBack(KratosMultiphysics.RIGID)
            node_flags.PushBack(KratosMultiphysics.NOT_FLUID)
            
            for fluid_part in fluid_body_model_parts:
                for rigid_part in rigid_body_model_parts:
                    transfer_process = TransferNodesProcess(fluid_part,rigid_part,node_flags)
                    transfer_process.Execute()
 
                            
            StopTimeMeasuring(clock_time,"1.rigid_body_model_parts  part.Nodes", True);

            
        #construct the computing model part:        
        domain_parts = []
        for i in range(self.sub_model_part_names.size()):
            domain_parts.append(self.main_model_part.GetSubModelPart(self.sub_model_part_names[i].GetString()))
        
        processes_parts = []
        for i in range(self.processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
       
        #construct a model part which contains the mesh to compute
        self.main_model_part.CreateSubModelPart(self.computing_model_part_name)
        fluid_computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        fluid_computing_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        fluid_computing_model_part.Properties  = self.main_model_part.Properties

        #set flag to identify the fluid model part
        fluid_computing_model_part.Set(KratosMultiphysics.FLUID)
        #set flag to identify the computing model part
        fluid_computing_model_part.Set(KratosMultiphysics.ACTIVE)
       
        for node in self.main_model_part.Nodes:
            #fluid_computing_model_part.AddNode(node,0)
            fluid_computing_model_part.Nodes.append(node)

        for part in domain_parts:
            #part.Set(KratosMultiphysics.FLUID)
            for elem in part.Elements:
                #fluid_computing_model_part.AddElement(elem,0)
                fluid_computing_model_part.Elements.append(elem)
            #for node in part.Nodes:
            #    fluid_computing_model_part.Nodes.append(node)
            
        for part in processes_parts:
            part.Set(KratosMultiphysics.BOUNDARY)
            for cond in part.Conditions:
                fluid_computing_model_part.AddCondition(cond,0)  
                #fluid_computing_model_part.Conditions.append(cond)  

        #delete body parts: (materials have to be already assigned)
        if( self.bodies_list == True ):
            for i in range(self.bodies_parts_list.size()):
               #get body parts
                body_parts_name_list = self.bodies_parts_list[i]["parts_list"]
                for j in range(body_parts_name_list.size()):
                    self.main_model_part.RemoveSubModelPart(body_parts_name_list[j].GetString())
                    print("::[Model_Prepare]::Body Part Removed:", body_parts_name_list[j].GetString())
                    
        #for part in domain_parts:
        #    self.main_model_part.RemoveSubModelPart(part.Name)

        print(self.main_model_part)
