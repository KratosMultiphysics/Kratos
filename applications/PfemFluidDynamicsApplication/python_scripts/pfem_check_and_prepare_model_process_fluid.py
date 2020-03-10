from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

import time as timer

def StartTimeMeasuring():
    """This function starts time calculation
    """
    # Measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    """This function ends time calculation
    """
    # Measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[PFEM_FLUID_MODEL]:: [ %.2f" % round(used_time,2),"s", process," ] ")


def Factory(settings, Model):
    """This function creates the process
    """
    if type(settings) != KratosMultiphysics.Parameters:
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class CheckAndPrepareModelProcess(KratosMultiphysics.Process):
    """The class of the CheckAndPrepareModelProcess class
    """
    def __init__(self, main_model_part, Parameters ):
        """The constructor of the CheckAndPrepareModelProcess-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        Keyword arguments:
        self -- It signifies an instance of a class.
        main_model_part -- The ModelPart to be used
        Parameters -- The settings for the process
        """
        self.main_model_part = main_model_part

        self.computing_model_part_name  = Parameters["computing_model_part_name"].GetString()
        if Parameters.Has("problem_domain_sub_model_part_list"):
            self.sub_model_part_names     = Parameters["problem_domain_sub_model_part_list"]
        if Parameters.Has("processes_sub_model_part_list"):
            self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

        self.bodies_parts_list = []
        self.bodies_list = False
        if Parameters.Has("bodies_list"):
            self.bodies_list = True
            self.bodies_parts_list = Parameters["bodies_list"]

        if Parameters.Has("material_import_settings"):
            self.material_import_settings = Parameters["material_import_settings"]

    def Execute(self):
        """This function executes the process
        """
        print("Execute check_and_prepare_model_process_fluid")

        #construct body model parts:
        fluid_body_model_parts = []
        solid_body_model_parts = []
        rigid_body_model_parts = []

        void_flags  = []

        if self.bodies_list == True:

            max_property_id = 0
            for prop in self.main_model_part.Properties:
                if prop.Id > max_property_id:
                    max_property_id = prop.Id

            for i in range(self.bodies_parts_list.size()):
                #create body model part
                body_model_part_name = self.bodies_parts_list[i]["body_name"].GetString()
                self.main_model_part.CreateSubModelPart(body_model_part_name)
                body_model_part = self.main_model_part.GetSubModelPart(body_model_part_name)

                print("[Model_Prepare]::Body Creation", body_model_part_name)
                body_model_part.ProcessInfo = self.main_model_part.ProcessInfo

                with open(self.material_import_settings["materials_filename"].GetString(), 'r') as parameter_file:
                    materials = KratosMultiphysics.Parameters(parameter_file.read())

                #build body from their parts
                body_parts_name_list = self.bodies_parts_list[i]["parts_list"]
                body_parts_list = []

                for j in range(body_parts_name_list.size()):

                    body_parts_list.append(self.main_model_part.GetSubModelPart(body_parts_name_list[j].GetString()))

                    # Add to the fluid and solid domains only the properties that belongs to their native SubModelParts
                    is_solid_or_fluid = False
                    for k in range(materials["properties"].size()):
                        if materials["properties"][k]["model_part_name"].GetString() == self.main_model_part.Name + "." + body_parts_name_list[j].GetString():
                            property_id = self.main_model_part.GetProperties()[materials["properties"][k]["properties_id"].GetInt()]
                            body_model_part.AddProperties(property_id)
                            is_solid_or_fluid = True

                    # Assign a new dummy property to the rigid boundaries to keep them separated in the GiD post process
                    if not is_solid_or_fluid:
                        max_property_id += 1
                        new_dummy_property = KratosMultiphysics.Properties(max_property_id)
                        self.main_model_part.GetSubModelPart(body_parts_name_list[j].GetString()).AddProperties(new_dummy_property)
                        KratosPfemFluid.SetDummyPropertyForRigidElementsProcess( \
                            self.main_model_part.GetSubModelPart(body_parts_name_list[j].GetString()), \
                            max_property_id).Execute()

                body_model_part_type = self.bodies_parts_list[i]["body_type"].GetString()

                for part in body_parts_list:
                    entity_type = "Nodes"
                    clock_time = StartTimeMeasuring()

                    if (body_model_part_type=="Fluid"):
                        assign_flags = [KratosMultiphysics.FLUID]
                        transfer_process = KratosDelaunay.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                        transfer_process.Execute()
                    elif (body_model_part_type=="Solid"):
                        assign_flags  = [KratosMultiphysics.SOLID]
                        transfer_process = KratosDelaunay.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                        transfer_process.Execute()
                    elif (body_model_part_type=="Rigid"):
                        assign_flags  = [KratosMultiphysics.RIGID,KratosMultiphysics.BOUNDARY]
                        transfer_process = KratosDelaunay.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                        transfer_process.Execute()

                    StopTimeMeasuring(clock_time,"1. for node in part.Nodes", True);

                    clock_time = StartTimeMeasuring()

                    entity_type = "Elements"
                    transfer_process = KratosDelaunay.TransferEntitiesProcess(body_model_part,part,entity_type)
                    transfer_process.Execute()

                    StopTimeMeasuring(clock_time,"1. part.Elements", True);

                    clock_time = StartTimeMeasuring()

                    entity_type = "Conditions"
                    transfer_process = KratosDelaunay.TransferEntitiesProcess(body_model_part,part,entity_type)

                    StopTimeMeasuring(clock_time,"1. part.Conditions", True);


                print(" Bodies Appended ")

                if body_model_part_type == "Fluid":
                    body_model_part.Set(KratosMultiphysics.FLUID)
                    fluid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
                if body_model_part_type == "Solid":
                    body_model_part.Set(KratosMultiphysics.FLUID)
                    body_model_part.Set(KratosMultiphysics.SOLID)
                    solid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
                if body_model_part_type == "Rigid":
                    body_model_part.Set(KratosMultiphysics.RIGID)
                    rigid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
            clock_time = StartTimeMeasuring()

            #add walls in fluid domains:
            transfer_flags = [KratosMultiphysics.RIGID, (KratosMultiphysics.FLUID).AsFalse()]

            for solid_part in solid_body_model_parts:
                set_solid_material_process=KratosPfemFluid.SetMaterialPropertiesToSolidNodes(solid_part)
                set_solid_material_process.Execute()

            entity_type = "Nodes"
            for fluid_part in fluid_body_model_parts:
                set_fluid_material_process=KratosPfemFluid.SetMaterialPropertiesToFluidNodes(fluid_part)
                set_fluid_material_process.Execute()

                for rigid_part in rigid_body_model_parts:
                    set_rigid_material_process=KratosPfemFluid.SetMaterialPropertiesFromFluidToRigidNodes(rigid_part,fluid_part)
                    set_rigid_material_process.Execute()
                    transfer_process = KratosDelaunay.TransferEntitiesProcess(fluid_part,rigid_part,entity_type,transfer_flags)
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

        entity_type = "Nodes"
        transfer_process = KratosDelaunay.TransferEntitiesProcess(fluid_computing_model_part,self.main_model_part,entity_type)
        transfer_process.Execute()

        for part in domain_parts:
            entity_type = "Elements"
            transfer_process = KratosDelaunay.TransferEntitiesProcess(fluid_computing_model_part,part,entity_type)
            transfer_process.Execute()

        for part in processes_parts:
            part.Set(KratosMultiphysics.BOUNDARY)
            entity_type = "Conditions"
            #condition flags as BOUNDARY or CONTACT are reserved to composite or contact conditions (do not set it here)
            transfer_process = KratosDelaunay.TransferEntitiesProcess(fluid_computing_model_part,part,entity_type)
            transfer_process.Execute()

        #delete body parts: (materials have to be already assigned)
        if self.bodies_list == True:
            for i in range(self.bodies_parts_list.size()):
               #get body parts
                body_parts_name_list = self.bodies_parts_list[i]["parts_list"]
                for j in range(body_parts_name_list.size()):
                    self.main_model_part.RemoveSubModelPart(body_parts_name_list[j].GetString())
                    print("::[Model_Prepare]::Body Part Removed:", body_parts_name_list[j].GetString())
        print(" Main Model Part", self.main_model_part )
