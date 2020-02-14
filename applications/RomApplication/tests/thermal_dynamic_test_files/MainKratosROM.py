from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.RomApplication as romapp

from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
from KratosMultiphysics.ConvectionDiffusionApplication.apply_thermal_face_process import ApplyThermalFaceProcess
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
from KratosMultiphysics.RomApplication import convection_diffusion_transient_rom_solver

import numpy as np
import json
import sys

class ConvDiffDynamicROM(ConvectionDiffusionAnalysis):

    def __init__(self,model,project_parameters):        
        super(ConvDiffDynamicROM,self).__init__(model,project_parameters)
        self.selected_time_step_solution_container = []

    def _CreateSolver(self):    
        return convection_diffusion_transient_rom_solver.ROMSolver(self.model, self.project_parameters["solver_settings"])
        
    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(ConvDiffDynamicROM,self).ModifyInitialGeometry()
        
        computing_model_part = self._solver.GetComputingModelPart()

        with open('NodalModes.json') as f:
            data=json.load(f)
            counter = 0
            rom_dofs= self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(1, rom_dofs) # 1, since its a scalar quantity (a vector is imported to each node)
                Counter=str(1.0*node.Id)
                for i in range(rom_dofs):
                    aux[0, i] = data[Counter][i] # ROM basis
                node.SetValue(romapp.ROM_BASIS, aux ) # Aux ID
                node.SetValue(romapp.AUX_ID, counter)
                counter+=1


    def ModifyInitialProperties(self):  
        self.processes = []
        #################################################################
        #                      HEAT FLUX PROPERTY                       #
        #################################################################    
        HeatFlux1 = 2000.0
        heatFluxSettings1 = KratosMultiphysics.Parameters("""
            {
                    "model_part_name": "ThermalModelPart.HeatFlux2D_Top_Wall",
                    "variable_name": "FACE_HEAT_FLUX",
                    "constrained": false,
                    "interval": [0.0,"End"]
            }
            """
            )
        heatFluxSettings1.AddEmptyValue("value").SetDouble(HeatFlux1)
        self.processes.append(AssignScalarVariableProcess(self.model, heatFluxSettings1))

        HeatFlux2 = 2000.0
        heatFluxSettings2 = KratosMultiphysics.Parameters("""
            {
                    "model_part_name": "ThermalModelPart.HeatFlux2D_Bottom_Wall",
                    "variable_name": "FACE_HEAT_FLUX",
                    "constrained": false,
                    "interval": [0.0,"End"]
            }
            """
            )
        heatFluxSettings2.AddEmptyValue("value").SetDouble(HeatFlux2)
        self.processes.append(AssignScalarVariableProcess(self.model, heatFluxSettings2))

        ################################################################
        #                     TEMPERATURE PROPERTY                     #
        ################################################################        
        ImposedTemperature = 303.15
        TemperatureSettings = KratosMultiphysics.Parameters("""
            {
                    "model_part_name": "ThermalModelPart.ImposedTemperature2D_Left_Wall",
                    "variable_name": "TEMPERATURE",
                    "constrained": true,
                    "interval": [0.0,"End"]
            }
            """
            )
        TemperatureSettings.AddEmptyValue("value").SetDouble(ImposedTemperature)
        self.processes.append(AssignScalarVariableProcess(self.model, TemperatureSettings))
        ################################################################

        # ################################################################
        # #                      RADIATION PROPERTY                      #
        # ################################################################
        
        RadiationSettings = KratosMultiphysics.Parameters("""
            {
                    "model_part_name": "ThermalModelPart.ThermalFace2D_Right_Wall",
                    "add_ambient_radiation": true,
                    "emissivity": 0.8,
                    "add_ambient_convection": true,
                    "convection_coefficient": 100.0,
                    "interval": [0.0,"End"]
            }
            """
            )

        AmbientTemperature = 283.15
        RadiationSettings.AddEmptyValue("ambient_temperature").SetDouble(AmbientTemperature)
        self.processes.append(ApplyThermalFaceProcess(self.model, RadiationSettings))

        ################################################################
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()
        ################################################################
    
    def FinalizeSolutionStep(self):
        super(ConvDiffDynamicROM,self).FinalizeSolutionStep()        
        ArrayOfTemperatures = []        
        if self.time==500 or self.time==1200 or self.time==2500 or self.time==3000 or self.time==3600:
            for node in self._solver.GetComputingModelPart().Nodes:       
                ArrayOfTemperatures.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
            self.selected_time_step_solution_container.append(ArrayOfTemperatures)  

    def EvaluateQuantityOfInterest(self):
       ################################################################################################
       # Functions evaluating the QoI of the problem: Numpy array with selected snapshots of solution #
       #                                    and nodal area                                            #
       ################################################################################################
        SnapshotMatrix = np.zeros((len(self.selected_time_step_solution_container[0]), len(self.selected_time_step_solution_container)))
        for i in range(len(self.selected_time_step_solution_container)):
            Snapshot_i= np.array(self.selected_time_step_solution_container[i])
            SnapshotMatrix[:,i] = Snapshot_i.transpose()
        return SnapshotMatrix

    def EvaluateQuantityOfInterest2(self):
        computing_model_part = self.model["ThermalModelPart"]
        dimension = self._GetSolver().settings["domain_size"].GetInt()
        area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(computing_model_part , dimension)
        area_calculator.Execute()        

        ArrayOfAreas = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        return ArrayOfAreas

############################################################################################################

if __name__ == "__main__":

    with open("ProjectParametersROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ConvDiffDynamicROM(model,parameters)
    simulation.Run()
