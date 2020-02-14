from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp

from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
from KratosMultiphysics.assign_vector_variable_process import AssignVectorVariableProcess
from KratosMultiphysics.assign_vector_by_direction_process import AssignVectorByDirectionProcess
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication import structural_mechanics_static_rom_solver

import numpy as np
import json
import sys

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

class StructStaticROM(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters):
        super(StructStaticROM,self).__init__(model,project_parameters)

    def _CreateSolver(self):        
        return structural_mechanics_static_rom_solver.ROMSolver(self.model, self.project_parameters["solver_settings"])
    
    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(StructStaticROM,self).ModifyInitialGeometry()
        
        computing_model_part = self._solver.GetComputingModelPart()

        with open('NodalModes.json') as f:
            data=json.load(f)
            counter = 0
            rom_dofs=self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            Dimensions = len(data['1.0'])  # dimensions are defined here (a vector or a matrix is imported to each node)
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(Dimensions, rom_dofs) 
                for j in range(Dimensions):
                    Counter=str(1.0*node.Id)
                    if Dimensions >1:
                        for i in range(rom_dofs):                         
                            aux[j,i] = data[Counter][j][i]
                    else:
                        for i in range(rom_dofs): 
                            aux[j,i] = data[Counter][i]
                node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                node.SetValue(romapp.AUX_ID, counter) # Aux ID
                counter+=1


    def ModifyInitialProperties(self):
        self.processes = []
        ############################################################################################
        #                      POINT LOAD CONDITION                                                #
        ############################################################################################
        
        PointLoad = 50.0
        PointLoadSettings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "Structure.PointLoad2D_Load_on_points_Auto1",
                "variable_name"   : "POINT_LOAD",
                "direction"       : [1.0,0.0,0.0],
                "interval"        : [0.0,"End"]
            }
            """
            )
        PointLoadSettings.AddEmptyValue("modulus").SetDouble(PointLoad)
        self.processes.append(AssignVectorByDirectionProcess(self.model, PointLoadSettings))



        ############################################################################################
        #                      LINE LOAD CONDITION                                                 #
        ############################################################################################

        LineLoad = 50.0
        LineLoadSettings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "Structure.LineLoad2D_InterfaceStructure",
                "variable_name"   : "LINE_LOAD",
                "direction"       : [1,0.0,0.0],
                "interval"        : [0.0,"End"]
            }
            """
            )
        LineLoadSettings.AddEmptyValue("modulus").SetDouble(LineLoad)
        self.processes.append(AssignVectorByDirectionProcess(self.model, LineLoadSettings))


        LineLoad2 = 50.0
        LineLoadSettings2 = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "Structure.LineLoad2D_Load_on_lines_Auto1",
                "variable_name"   : "LINE_LOAD",
                "direction"       : [1.0,0.0,0.0],
                "interval"        : [0.0,"End"]
            }
            """
            )
        LineLoadSettings2.AddEmptyValue("modulus").SetDouble(LineLoad2)
        self.processes.append(AssignVectorByDirectionProcess(self.model, LineLoadSettings2))


        ############################################################################################
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()
        ############################################################################################

    def EvaluateQuantityOfInterest(self):
       ##############################################################################################
       # Functions evaluating the QoI of the problem: Array of displacement at every node on the mesh #
       #                    and nodal area                                                           #
       ##############################################################################################
        ArrayOfDisplacements = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:            
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0))
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0))
        return ArrayOfDisplacements

    def EvaluateQuantityOfInterest2(self):
        computing_model_part = self._solver.GetComputingModelPart()
        dimension = self._GetSolver().settings["domain_size"].GetInt()
        area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(computing_model_part , dimension)
        area_calculator.Execute()        

        ArrayOfAreas = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        return ArrayOfAreas    


if __name__ == "__main__":
    with open("ProjectParametersROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()      
    Simulation = StructStaticROM(model,parameters)
    Simulation.Run()
