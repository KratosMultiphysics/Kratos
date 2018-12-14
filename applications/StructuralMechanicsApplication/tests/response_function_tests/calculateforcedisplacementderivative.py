# Import Kratos core and apps
import os, shutil
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import structural_response_function_factory
import structural_mechanics_analysis
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosHDF5Application
import h5py
import numpy as np
import matplotlib.pyplot as plt
import time as timer
from gid_output_process import GiDOutputProcess
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication

with open("nonlinear_shell_cantilever_surface_load.json",'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters( parameter_file.read())

model = KratosMultiphysics.Model()

model_part_primal = model.CreateModelPart("Structure" , 2)
primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, parameters)
primal_analysis.Initialize()

startTime = timer.time()
if not primal_analysis.time < primal_analysis.end_time:
    primal_analysis.end_time += 1

primal_analysis.RunSolutionLoop()

gid_output = GiDOutputProcess(model_part_primal,
                             "gid_output",
                            KratosMultiphysics.Parameters("""
                                {
                                    "result_file_configuration" : {
                                        "gidpost_flags": {
                                            "GiDPostMode": "GiD_PostBinary",
                                            "WriteDeformedMeshFlag": "WriteUndeformed",
                                            "WriteConditionsFlag": "WriteConditions",
                                            "MultiFileFlag": "SingleFile"
                                        },
                                        "file_label"          : "step",
                                        "output_control_type" : "step",
                                         "output_frequency"    : 1,
                                        "nodal_results"       : ["DISPLACEMENT" , "ROTATION"]
                                    }
                                }
                                """)
                            )

gid_output.ExecuteInitialize()
gid_output.ExecuteBeforeSolutionLoop()
gid_output.ExecuteInitializeSolutionStep()
gid_output.PrintOutput()
gid_output.ExecuteFinalizeSolutionStep()
gid_output.ExecuteFinalize()


LHS_1 = KratosMultiphysics.Matrix(6,6)
RHS_1 = KratosMultiphysics.Vector(12)

RHS_unperturbed = []
for condition in model_part_primal.Conditions:
    condition.CalculateLocalSystem(LHS_1,RHS_1,model_part_primal.ProcessInfo)
    RHS_unperturbed.append(RHS_1)
    print("RHS_unperturbed" ,  RHS_1)

i = 1
perturbation_vector_X = []
perturbed_RHS_X = []
perturbed_RHS_Y = []
perturbed_RHS_Z = []
RHS_difference_X = []
RHS_difference_Y = []
RHS_difference_Z = []

#calculation of the difference in RHS due to displacement perturbation
for node in model_part_primal.Nodes:
    RHS_difference_X_node = []
    RHS_difference_Y_node = []
    RHS_difference_Z_node = []
    LHS = KratosMultiphysics.Matrix(6,6)
    RHS = KratosMultiphysics.Vector(12)
    displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0)
    #rotation = node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
    perturbation_disp = 0.001 * displacement
    node.X += perturbation_disp[0]
    for condition in model_part_primal.Conditions:
        condition.CalculateLocalSystem(LHS,RHS,model_part_primal.ProcessInfo)
        perturbed_RHS_X.append(RHS)
        RHS_difference_X_node.append(RHS - RHS_unperturbed[condition.Id - 1])
    node.X -= perturbation_disp[0]
    print("RHS_difference_X_node" , "node Id" , node.Id , ":" , RHS_difference_X_node[9])

    LHS_2 = KratosMultiphysics.Matrix(6,6)
    RHS_2 = KratosMultiphysics.Vector(12)
    node.Y += perturbation_disp[1]
    for condition in model_part_primal.Conditions:
        condition.CalculateLocalSystem(LHS_2,RHS_2,model_part_primal.ProcessInfo)
        perturbed_RHS_Y.append(RHS_2)
        RHS_difference_Y_node.append(RHS_2 - RHS_unperturbed[condition.Id - 1])
    node.Y -= perturbation_disp[1]
    print("RHS_difference_Y_node" , "node Id" , node.Id , ":" , RHS_difference_Y_node[9])

    LHS_3 = KratosMultiphysics.Matrix(6,6)
    RHS_3 = KratosMultiphysics.Vector(12)
    node.Z += perturbation_disp[2]
    for condition in model_part_primal.Conditions:
        condition.CalculateLocalSystem(LHS_3,RHS_3,model_part_primal.ProcessInfo)
        perturbed_RHS_Z.append(RHS_3)
        RHS_difference_Z_node.append(RHS_3 - RHS_unperturbed[condition.Id - 1])
    node.Z -= perturbation_disp[2]
    print("RHS_difference_Z_node" , "node Id" , node.Id , ":" , RHS_difference_Z_node[9])

    RHS_difference_X.append(RHS_difference_X_node)
    RHS_difference_Y.append(RHS_difference_Y_node)
    RHS_difference_Z.append(RHS_difference_Z_node)

    i+=1



