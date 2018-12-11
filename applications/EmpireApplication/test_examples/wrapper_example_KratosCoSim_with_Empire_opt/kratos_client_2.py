from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
# EMPIRE
import KratosMultiphysics.EmpireApplication as KratosEmpire
from ctypes import cdll
import os
import ctypes as ctp

def SetNodalValues(counter): # Somehow Modify the Nodal Values
    for node in model_part.Nodes:
        value = sum(node.GetSolutionStepValue(VELOCITY)) + counter
        node.SetSolutionStepValue(PRESSURE, value)

print("This is kratos_client_2")

model_part = ModelPart("MyModelPart")
model_part.AddNodalSolutionStepVariable(VELOCITY)
model_part.AddNodalSolutionStepVariable(PRESSURE)

print("Starting to initialize Empire")
import empire_wrapper
print("Import Successfull")
empire = empire_wrapper.EmpireWrapper()
print("Wrapper Created")
empire.Connect("kratos_client_2.xml")

empire.ReceiveMesh("myMesh2", model_part)

for i in range(10):
    is_converged = 0
    k = 0
    while not is_converged:
        k += 1
        empire.ReceiveDataField("myMesh2", "velocity", [VELOCITY, VELOCITY]) # receiving an (Empire) "doubleVector"
        SetNodalValues(i)
        empire.SendDataField("myMesh2", "pressure", PRESSURE)

        is_converged = empire.ReceiveConvergenceSignal()
        print("kratos_client_2; i =", i, "| k =", k, "| is_converged:", is_converged)

empire.Disconnect()