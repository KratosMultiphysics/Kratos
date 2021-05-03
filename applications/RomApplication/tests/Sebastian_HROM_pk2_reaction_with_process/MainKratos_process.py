import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
from matplotlib import pyplot as plt
import numpy as np
import json
import time

import pdb

class RunHROM(StructuralMechanicsAnalysisROM):

    def ModifyInitialGeometry(self):
        """Here is the place where the HROM_WEIGHTS are assigned to the selected elements and conditions"""
        super().ModifyInitialGeometry()
        computing_model_part = self._solver.GetComputingModelPart()
        ## Adding the weights to the corresponding elements
        with open('ElementsAndWeights.json') as f:
            HR_data = json.load(f)
            for key in HR_data["Elements"].keys():
                computing_model_part.GetElement(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Elements"][key])
            for key in HR_data["Conditions"].keys():
                computing_model_part.GetCondition(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Conditions"][key])

if __name__ == "__main__":
    ##############################################################################################
    #                                           TRAIN ROM                                        #
    ##############################################################################################
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model,parameters)
    simulation.Run()

    ##############################################################################################
    #                                          TRAIN HROM                                        #
    ##############################################################################################
    start = time.perf_counter()
    with open("ProjectParameters_ROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysisROM(model,parameters,"EmpiricalCubature")
    simulation.Run()
    finish = time.perf_counter()
    print("Finished in ",round(finish-start,2)," second(s)")

    ##############################################################################################
    #                                           RUN HROM                                         #
    ##############################################################################################
    start = time.perf_counter()
    with open("ProjectParameters_HROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = RunHROM(model,parameters)
    simulation.Run()
    finish = time.perf_counter()
    print("Finished in ",round(finish-start,2)," second(s)")