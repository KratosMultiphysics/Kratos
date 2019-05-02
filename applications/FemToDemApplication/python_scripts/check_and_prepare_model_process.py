from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
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

        void_flags = []

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
                    entity_type = "Nodes"
                    if (body_model_part_type=="Fluid"):
                        assign_flags = [KratosMultiphysics.FLUID]
                        transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                        transfer_process.Execute()
                    elif (body_model_part_type=="Solid"):
                        assign_flags = [KratosMultiphysics.SOLID]
                        transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                        transfer_process.Execute()
                    elif (body_model_part_type=="Rigid"):
                        assign_flags = [KratosMultiphysics.RIGID,KratosMultiphysics.BOUNDARY]
                        transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                        transfer_process.Execute()

                    entity_type = "Elements"
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type)
                    transfer_process.Execute()
                    entity_type = "Conditions"
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type)
                    transfer_process.Execute()

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
            transfer_flags = [KratosMultiphysics.RIGID,KratosMultiphysics.NOT_FLUID]

            entity_type = "Nodes"
            for fluid_part in fluid_body_model_parts:
                for rigid_part in rigid_body_model_parts:
                    transfer_process = KratosSolid.TransferEntitiesProcess(fluid_part,rigid_part,entity_type,transfer_flags)
                    transfer_process.Execute()


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


        entity_type = "Nodes"
        transfer_process = KratosSolid.TransferEntitiesProcess(computing_model_part,self.main_model_part,entity_type)
        transfer_process.Execute()

        for part in domain_parts:
            entity_type = "Elements"
            transfer_process = KratosSolid.TransferEntitiesProcess(computing_model_part,part,entity_type)
            transfer_process.Execute()

        for part in processes_parts:
            part.Set(KratosMultiphysics.BOUNDARY)
            entity_type = "Conditions"
            #condition flags as BOUNDARY or CONTACT are reserved to composite or contact conditions (do not set it here)
            transfer_process = KratosSolid.TransferEntitiesProcess(computing_model_part,part,entity_type)
            transfer_process.Execute()

        #delete body parts: (materials have to be already assigned)
        if( self.bodies_list == True ):
            for i in range(self.bodies_parts_list.size()):
               #get body parts
                body_parts_name_list = self.bodies_parts_list[i]["parts_list"]
                for j in range(body_parts_name_list.size()):
                    self.main_model_part.RemoveSubModelPart(body_parts_name_list[j].GetString())
                    print("::[Model_Prepare]::Body Part Removed:", body_parts_name_list[j].GetString())
