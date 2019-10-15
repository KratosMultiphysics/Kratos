from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.PfemFluidDynamicsApplication.pfem_check_and_prepare_model_process_fluid as pfem_check_and_prepare_model_process_fluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import time as timer

def Wait():
    input("CheckAndPrepareModelProcess -> Press Something")

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
        print("::[PFEM_FLUID_MODEL]:: [ %.2f" % round(used_time, 2), "s", process, " ] ")

#============================================================================================================================
class CheckAndPrepareModelProcessForCoupling(pfem_check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess):
#============================================================================================================================
    """
    The class of the CheckAndPrepareModelProcess class
    """
#============================================================================================================================
    def __init__(self, main_model_part, FEM_model_part, Parameters ):
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
        self.FEM_model_part  = FEM_model_part

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

#============================================================================================================================
    def Execute(self):
        """This function executes the process
        """
        self.AddAndReorderFEMDEMBoundary()
        super(CheckAndPrepareModelProcessForCoupling, self).Execute()

#============================================================================================================================
    def AddAndReorderFEMDEMBoundary(self):
        
        self.FEM_model_part.CreateSubModelPart("SkinDEMModelPart")
        skin_params = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 1
        }""")

        skin_detection_process = KratosMultiphysics.SkinDetectionProcess2D(self.FEM_model_part, skin_params)
        skin_detection_process.Execute()

        max_id = 0
        for node in self.main_model_part.Nodes:
            if node.Id > max_id:
                max_id = node.Id
        self.main_model_part.CreateSubModelPart("FEMDEM_boundary")
        femdem_model_part = self.main_model_part.GetSubModelPart("FEMDEM_boundary")

        # Reorder nodes Id
        for node in self.FEM_model_part.Nodes:
            node.Id = node.Id + max_id
            node.SetValue(KratosPfemFluid.NO_MESH, True)
            femdem_model_part.AddNode(node, 0)

        for node in self.FEM_model_part.GetSubModelPart("SkinDEMModelPart").Nodes:
            node.SetValue(KratosPfemFluid.NO_MESH, False)
#============================================================================================================================