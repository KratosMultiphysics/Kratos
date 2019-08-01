from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import h5py
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics import python_solver
import json

from convection_diffusion_analysis import ConvectionDiffusionAnalysis

import sys
import time

class ConvectionDiffusionAnalysisWithFlush(ConvectionDiffusionAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(ConvectionDiffusionAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def _CreateSolver(self):
        print(self.parallel_type)
        if self.parallel_type == "OpenMP":
            import rom_solver
            return rom_solver.ROMSolver(self.model, self.project_parameters["solver_settings"])
        else:
            import trilinos_static_solver
            return trilinos_static_solver.TrilinosStaticThermalSolver(self.model, self.project_parameters["solver_settings"])

    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(ConvectionDiffusionAnalysisWithFlush,self).ModifyInitialGeometry()
        
        computing_model_part = self.model["ThermalModelPart"]
        
        with open('NewProjectParameters.json') as f:
            data=json.load(f)
            counter = 0
            rom_dofs=self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(1, rom_dofs)
                Counter=str(1.0*node.Id)
                for i in range(rom_dofs):                
                    aux[0, i] = data[Counter][i]
                """ROM_BASIS"""    
                node.SetValue(KratosMultiphysics.ROM_BASIS, aux )
                """AUX_ID"""
                node.SetValue(KratosMultiphysics.AUX_ID, counter)
                counter+=1

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

    def ModifyInitialProperties(self):       
     

        # ################################################################
        # #                      HEAT FLUX PROPERTY                      #
        ################################################################  
        model_part = self.model.GetModelPart("ThermalModelPart.HeatFlux2D_Top_Wall")
        mesh = model_part.GetMesh(0)
        is_fixed = True
        value = -700.0
        variable = KratosMultiphysics.KratosGlobals.GetVariable("FACE_HEAT_FLUX")
        variable_utils = KratosMultiphysics.VariableUtils()
        variable_utils.ApplyFixity(variable, is_fixed, mesh.Nodes)
        variable_utils.SetScalarVar(variable, value, mesh.Nodes)
        
        
        model_part = self.model.GetModelPart("ThermalModelPart.HeatFlux2D_Bottom_Wall")
        mesh = model_part.GetMesh(0)
        is_fixed = True
        value = 700.0
        variable = KratosMultiphysics.KratosGlobals.GetVariable("FACE_HEAT_FLUX") ##Redundant
        variable_utils = KratosMultiphysics.VariableUtils()  ##Redundant
        variable_utils.ApplyFixity(variable, is_fixed, mesh.Nodes)
        variable_utils.SetScalarVar(variable, value, mesh.Nodes)    

           
                
        ################################################################
        
        ################################################################
        #                     TEMPERATURE PROPERTY                     #
        ################################################################ 
        model_part = self.model.GetModelPart("ThermalModelPart.ImposedTemperature2D_Left_Wall")
        mesh = model_part.GetMesh(0)
        is_fixed = True
        T = 90.0  

        # ## To use with function-dependent constraint
        # variable = KratosMultiphysics.KratosGlobals.GetVariable("TEMPERATURE")
        # aux_function = KratosMultiphysics.PythonGenericFunctionUtility( "(x**2 + y**2) * sin(0.5*t) " )
        # cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility( mesh.Nodes, aux_function )
        # variable_utils = KratosMultiphysics.VariableUtils()
        # variable_utils.ApplyFixity( variable, is_fixed, mesh.Nodes)
        # cpp_apply_function_utility.ApplyFunction( variable, time )


        # To use with constant constraint
        value = T
        variable = KratosMultiphysics.KratosGlobals.GetVariable("TEMPERATURE")
        variable_utils = KratosMultiphysics.VariableUtils()
        variable_utils.ApplyFixity(variable, is_fixed, mesh.Nodes)
        variable_utils.SetScalarVar(variable, value, mesh.Nodes)  
        
        
        ################################################################  

        
        # ################################################################
        # #                      RADIATION PROPERTY                      #
        # ################################################################        
        max_prop_id = -1
        model_part_names = self.model.GetModelPartNames()
        model_part = self.model.GetModelPart("ThermalModelPart")
        for model_part_name in model_part_names:
            for prop in self.model.GetModelPart(model_part_name).Properties:
                if prop.Id > max_prop_id:
                    max_prop_id = prop.Id
        model_part.GetCommunicator().MaxAll(max_prop_id)

        # Create a new property with the user defined interface parameters
        thermal_interface_prop = KratosMultiphysics.Properties(max_prop_id + 1)
        emissivity = 0.8
        convection_coefficient = 100.0
        ambient_temperature = 12.0
        thermal_interface_prop.SetValue(KratosMultiphysics.EMISSIVITY, emissivity)
        thermal_interface_prop.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, convection_coefficient)
        thermal_interface_prop.SetValue(KratosMultiphysics.AMBIENT_TEMPERATURE, ambient_temperature)

        # Set the new property in the thermal face model part
        model_part = self.model.GetModelPart("ThermalModelPart.ThermalFace2D_Right_Wall")
        model_part.AddProperties(thermal_interface_prop)
        for condition in model_part.Conditions:
            condition.Properties = thermal_interface_prop
            
            
        #####################################################################
      
       ##############################################################################################
       # function evaluating the QoI of the problem: Array of temperature at every node on the mesh #
       ##############################################################################################
    def EvaluateQuantityOfInterest(self):
        ArrayOfTemperatures = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfTemperatures.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        return ArrayOfTemperatures


        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now
                
        ##############################################################################################                

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
        
    model = KratosMultiphysics.Model()
    simulation = ConvectionDiffusionAnalysisWithFlush(model,parameters)
    simulation.Run()
    QoI = simulation.EvaluateQuantityOfInterest()
    print( "------------------------------------------" )
    print("Temperature list saved in hdf5 format")

    with h5py.File('./ROM.h5','w') as hdf:
        hdf.create_dataset('Temperature',data=QoI)
