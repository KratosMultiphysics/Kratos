import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
from matplotlib import pyplot as plt
import numpy as np
import json
import time

import pdb

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
    with open("ProjectParameters_ROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysisROM(model,parameters,"EmpiricalCubature")
    simulation.Run()
