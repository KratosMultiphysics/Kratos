from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
# EMPIRE
import KratosMultiphysics.EmpireApplication as KratosEmpire
from ctypes import cdll
import os
import ctypes as ctp
import numpy as np

def SetNodalValues(counter): # Somehow Modify the Nodal Values
    for node in model_part.Nodes:
        value = node.GetSolutionStepValue(TEMPERATURE) - counter**2
        node.SetSolutionStepValue(DISPLACEMENT, Vector([value, value+counter, value+counter*2]))

print("This is kratos_client_1")

model_part = ModelPart("MyModelPart")
model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(TEMPERATURE)

nodes = np.array([0.0,0.0,0.0, 1.0,0.0,0.0, 2.0,0.0,0.0, 3.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, 2.0,1.0,0.0, 3.0,1.0,0.0])
elements = np.array([1,2,6,5, 2,3,7,6, 3,4,8,7])

# Creating the mesh in the model part
print("Mesh 1 ::::: Length of Nodes = ",int(len(nodes)/3))
print("Mesh 1 ::::: Length of Elements = ",int(len(elements)/4))
print("Mesh 1 ::::: Maximum entry in elements vector = ", np.max(elements))
for i in range(0,int(len(nodes)/3)):
    model_part.CreateNewNode(i+1, nodes[3*i+0], nodes[3*i+1], nodes[3*i+2])
print("Mesh 1 ::::: Finished adding nodes to model part !!")

# Creating elements with the above nodes
model_part.AddProperties(Properties(1))
print("Mesh 1 ::::: Finished adding Properties to model part !!")
prp = model_part.GetProperties()[1]
print("Mesh 1 ::::: Finished extracting properties from the model part !!")

for i in range(0,int(len(elements)/4)):
    sys.stdout.flush()
    model_part.CreateNewElement("Element2D4N", i+1,[int(elements[4*i+0]),int(elements[4*i+1]),int(elements[4*i+2]),int(elements[4*i+3])], prp)
print("Mesh 1 ::::: Finished Adding Elements to model part !!")

print("Starting to initialize Empire")
import empire_wrapper
print("Import Successfull")
empire = empire_wrapper.EmpireWrapper(echo_level=2)
print("Wrapper Created")
empire.Connect("kratos_client_1.xml")

empire.SendMesh("myMesh1", model_part)

for i in range(10):
    empire.SendDataField("myMesh1", "displacement", [DISPLACEMENT, DISPLACEMENT]) # sending an (Empire) "doubleVector"
    empire.ReceiveDataField("myMesh1", "temperature", TEMPERATURE)

    SetNodalValues(i)

empire.Disconnect()