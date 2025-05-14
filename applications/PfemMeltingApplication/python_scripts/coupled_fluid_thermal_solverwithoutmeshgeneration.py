from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
import KratosMultiphysics.MeshingApplication as MeshApp
import KratosMultiphysics.PfemMeltingApplication as PfemM

import time as timer

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):

    return PfemCoupledFluidThermalSolver(main_model_part, custom_settings)

import KratosMultiphysics.PfemMeltingApplication.coupled_fluid_thermal_pfem2solver
BaseClass = KratosMultiphysics.PfemMeltingApplication.coupled_fluid_thermal_pfem2solver.PfemCoupledFluidThermalSolver

class PfemCoupledFluidThermalSolver(BaseClass):
    def GetSuffix(self):
        return "_no_remesh"

    def PrepareModelPart(self):
        self.fluid_solver.PrepareModelPart()
        self.thermal_solver.PrepareModelPart()

        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 0.0)
        parametersf=self.settings["fluid_solver_settings"]

        self.cleaning_submodelparts()
        self.assign_nodally_properties();
        self.ReMesh()

    def cleaning_submodelparts(self):
        parametersf=self.settings["fluid_solver_settings"]
        parameters=self.settings["thermal_solver_settings"]
        self.skin_parts_list = []
        if parametersf.Has("skin_parts"):
            self.skin_parts_list = parametersf["skin_parts"]
        for i in range(self.skin_parts_list.size()):
            body_model_part_name=self.skin_parts_list[i].GetString()
            body_model_part_name=self.fluid_solver.main_model_part.GetSubModelPart(body_model_part_name)
            body_model_part_name.Conditions.clear()
            body_model_part_name.Elements.clear()
            body_model_part_name.Nodes.clear()

        self.bodies_parts_list = []
        if parameters.Has("processes_sub_model_part_list"):
            self.bodies_parts_list = parameters["processes_sub_model_part_list"]

    def filling_submodelparts(self):
        fluid_computational_model_part=self.fluid_solver.main_model_part.GetSubModelPart("fluid_computational_model_part")
        fluid_computational_model_part.Conditions.clear()
        fluid_computational_model_part.Elements.clear()
        fluid_computational_model_part.Nodes.clear()

        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(fluid_computational_model_part, self.fluid_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
        transfer_process.Execute()

        fluid_computational_model_part.ProcessInfo = self.fluid_solver.main_model_part.ProcessInfo
        fluid_computational_model_part.Properties  = self.fluid_solver.main_model_part.Properties
        if not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain"):
            self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")

        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part, "EulerianConvDiffLumped2D", "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part,"EulerianConvDiffLumped3D","ThermalFace3D3N")

    def filling_submodelparts_aux(self):
        if self.domain_size == 2:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part, "EulerianConvDiffLumped2D", "ThermalFace2D2N")
        else:
            self.modeler.GenerateModelPart(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part,"EulerianConvDiffLumped3D","ThermalFace3D3N")

    def ReMesh(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,self.mesh_element_size.GetDouble());

        for node in (self.fluid_solver.main_model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)

        self.fluid_solver.main_model_part.Conditions.clear()
        self.fluid_solver.main_model_part.Elements.clear()

        (self.Mesher).ReGenerateMesh("LagrangianFluidVMS3D","ThermalFace3D3N", self.fluid_solver.main_model_part, self.node_erase_process, True, False, 1.4, 0.1)  #1.8

        #LagrangianFluidVMS3D
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        neighbor_search = KratosMultiphysics.FindGlobalNodalNeighboursProcess(kratos_comm,self.fluid_solver.main_model_part)

        neighbor_search.Execute()
        neighbor_elements_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(kratos_comm,self.fluid_solver.main_model_part)

        neighbor_elements_search.Execute()
        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()

        (self.PfemM_apply_bc_process).Execute();
        self.node_erase_process.Execute()

        skin_params = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 1
        }""")

        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.fluid_solver.main_model_part.Nodes)
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()
        normal_calculation_utils.CalculateOnSimplex(self.fluid_solver.main_model_part.Conditions, self.domain_size)
        for node in self.fluid_solver.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY)== 1.0):
                solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
                solution_normal /= math.sqrt(solution_normal[0]**2+solution_normal[1]**2+solution_normal[2]**2)
                node.SetSolutionStepValue(KratosMultiphysics.NORMAL, solution_normal)

    def AuxReMesh(self):
        for node in (self.fluid_solver.main_model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)

        self.fluid_solver.main_model_part.Conditions.clear()

        neighbor_search = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.fluid_solver.main_model_part)
        neighbor_search.Execute()

        (self.PfemM_apply_bc_process).Execute();

        skin_params = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 0
        }""")

        skin_detection_process = KratosMultiphysics.SkinDetectionProcess3D(self.fluid_solver.main_model_part, skin_params)
        skin_detection_process.Execute()

        neighbor_condition_search = KratosMultiphysics.FindConditionsNeighboursProcess(self.fluid_solver.main_model_part,3, 20)
        neighbor_condition_search.Execute()

        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY,0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_FREE_SURFACE,0, 0.0)

        for node in self.fluid_solver.main_model_part.GetSubModelPart("SkinDEMModelPart").Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY,0, 1.0)

        (self.PfemM_apply_bc_process).Execute();

        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.fluid_solver.main_model_part.Nodes)
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()

        normal_calculation_utils.CalculateOnSimplex(self.fluid_solver.main_model_part.Conditions, self.domain_size)
        for node in self.fluid_solver.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.IS_BOUNDARY)== 1.0):
                solution_normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
                solution_normal /= math.sqrt(solution_normal[0]**2+solution_normal[1]**2+solution_normal[2]**2)
                node.SetSolutionStepValue(KratosMultiphysics.NORMAL, solution_normal)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.FACE_HEAT_FLUX, 0.0)

        return new_time

    def InitializeSolutionStep(self):
        self.step=1
        if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):

            for node in self.fluid_solver.main_model_part.Nodes:
                if(node.IsFixed(KratosMultiphysics.VELOCITY_X)==True):
                    node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0, 1.0) #NODES NOT
                    node.Fix(KratosMultiphysics.TEMPERATURE)

        t1 = timer.time()
        self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)

        t2= timer.time()
        self.streamlineintegration=self.streamlineintegration + t2 - t1

        self.inverted=False

        self.inverted=self.Streamline.CheckInvertElement(self.fluid_solver.main_model_part,self.domain_size,self.mesh_element_size.GetDouble())

        self.AuxReMesh()

        t3=timer.time()
        self.meshingprocedure=self.meshingprocedure + t3 - t2

        #Volume evaluation
        volume=self.Streamline.CalculateVolume(self.fluid_solver.main_model_part,3)
        step=self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        self.outputfile4.write(str(step)+" "+ str(volume) +"\n")

        t4=timer.time()
        if(self.step==self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]):

            self.filling_submodelparts()
        else:
            self.filling_submodelparts_aux()

        t5=timer.time()
        self.fillingsubmodelparts=self.fillingsubmodelparts + t5 - t4


        self.fluid_solver.InitializeSolutionStep()

        self.thermal_solver.InitializeSolutionStep()

        t6=timer.time()
        self.initializeSolutionStep=self.initializeSolutionStep + t6 - t5

