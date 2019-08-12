from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.RomApplication as romapp

from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
from KratosMultiphysics.ConvectionDiffusionApplication.apply_thermal_face_process import ApplyThermalFaceProcess
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

import h5py
import json
import sys
import time

class ConvectionDiffusionAnalysisWithFlush(ConvectionDiffusionAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(ConvectionDiffusionAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def _CreateSolver(self):
        print(sys.path, file=sys.stderr)
        print(self.parallel_type)
        if self.parallel_type == "OpenMP":
            from KratosMultiphysics.RomApplication import ConvectionDiffusionROMsolver
            return ConvectionDiffusionROMsolver.ROMSolver(self.model, self.project_parameters["solver_settings"])
        else:
            import trilinos_static_solver
            return trilinos_static_solver.TrilinosStaticThermalSolver(self.model, self.project_parameters["solver_settings"])

    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(ConvectionDiffusionAnalysisWithFlush,self).ModifyInitialGeometry()
        
        computing_model_part = self.model["ThermalModelPart"]
        
        with open('NodalModes.json') as f:
            data=json.load(f)
            counter = 0
            rom_dofs=self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(1, rom_dofs) # 1, since its a scalar quantity (a vector is imported to each node)
                Counter=str(1.0*node.Id)
                for i in range(rom_dofs):
                    aux[0, i] = data[Counter][i]
                """ROM_BASIS"""
                node.SetValue(romapp.ROM_BASIS, aux )
                """AUX_ID"""
                node.SetValue(romapp.AUX_ID, counter)
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
        self.processes = []

        HeatFlux1 = -700.0
        heatFluxSettings1 = KratosMultiphysics.Parameters("""
            {
                "model_part_name"             : "ThermalModelPart.HeatFlux2D_Top_Wall",
                "variable_name"               : "FACE_HEAT_FLUX",
                "constrained"                 : true,
                "interval"                    : [0.0,"End"]
            }
            """
            )
        heatFluxSettings1.AddEmptyValue("value").SetDouble(HeatFlux1)
        self.processes.append(AssignScalarVariableProcess(self.model, heatFluxSettings1))

        
        HeatFlux2 = 700.0
        heatFluxSettings2 = KratosMultiphysics.Parameters("""
            {
                "model_part_name"             : "ThermalModelPart.HeatFlux2D_Bottom_Wall",
                "variable_name"               : "FACE_HEAT_FLUX",
                "constrained"                 : true,
                "interval"                    : [0.0,"End"]
            }
            """
            )
        heatFluxSettings2.AddEmptyValue("value").SetDouble(HeatFlux2)
        self.processes.append(AssignScalarVariableProcess(self.model, heatFluxSettings2))

        ################################################################
        #                     TEMPERATURE PROPERTY                     #
        ################################################################
        #time = 1.0
        ImposedTemperature = 90
        TemperatureSettings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"             : "ThermalModelPart.ImposedTemperature2D_Left_Wall",
                "variable_name"               : "TEMPERATURE",
                "constrained"                 : true,
                "interval"                    : [0.0,"End"]
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
                "convection_coefficient": 100.0
            }
            """
            )

        AmbientTemperature = 12.0
        RadiationSettings.AddEmptyValue("ambient_temperature").SetDouble(AmbientTemperature)
        self.processes.append(ApplyThermalFaceProcess(self.model, RadiationSettings))
        ################################################################
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()
        ################################################################

    def EvaluateQuantityOfInterest(self):
       ##############################################################################################
       # Functions evaluating the QoI of the problem: Array of temperature at every node on the mesh #
       #                    and nodal area                                                           #
       ##############################################################################################
        ArrayOfTemperatures = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfTemperatures.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        return ArrayOfTemperatures

    def EvaluateQuantityOfInterest2(self):
        ArrayOfAreas = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        return ArrayOfAreas

    def FinalizeSolutionStep(self):
        super(ConvectionDiffusionAnalysisWithFlush,self).FinalizeSolutionStep()
        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

        ##############################################################################################

if __name__ == "__main__":

    with open("ProjectParametersROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ConvectionDiffusionAnalysisWithFlush(model,parameters)
    simulation.Run()
    
    QoI2 = simulation.EvaluateQuantityOfInterest2()
    print(QoI2)
    print( "------------------------------------------" )
    print("Area list saved in hdf5 format")

    with h5py.File('./AREA.h5','w') as hdf:
        hdf.create_dataset('Area',data=QoI2)
    print( "------------------------------------------" )

    QoI = simulation.EvaluateQuantityOfInterest()
    print(QoI)
    print( "------------------------------------------" )
    print("Temperature list saved in hdf5 format")

    with h5py.File('./FINE.h5','w') as hdf:
        hdf.create_dataset('Temperature',data=QoI)
