# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.MappingApplication import *
from KratosMultiphysics.StructuralApplication import *
import RBFMapper as RBFMapper

import numpy as np


# defining a model part for the fluid
mesh1_model_part = ModelPart("mesh1_coarse")
mesh2_model_part = ModelPart("mesh2_fine")

mesh1_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
mesh1_model_part.AddNodalSolutionStepVariable(PATCH_INDEX)

mesh2_model_part.AddNodalSolutionStepVariable(VELOCITY)
mesh2_model_part.AddNodalSolutionStepVariable(PATCH_INDEX)

# Mesh 1  Definition
nodes = np.array([0.0,0.0,0.0, 1.0,0.0,0.0, 2.0,0.0,0.0, 3.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, 2.0,1.0,0.0, 3.0,1.0,0.0])
elements = np.array([1,2,6,5, 2,3,7,6, 3,4,8,7])

# Creating the mesh in the model part
print("Mesh 1 ::::: Length of Nodes = ",int(len(nodes)/3))
print("Mesh 1 ::::: Length of Elements = ",int(len(elements)/4))
for i in range(0,int(len(nodes)/3)):
    mesh1_model_part.CreateNewNode(i+1, nodes[3*i+0], nodes[3*i+1], nodes[3*i+2])

# Creating elements with the above nodes
mesh1_model_part.AddProperties(Properties(1))
prp = mesh1_model_part.GetProperties()[1]

for node in mesh1_model_part.Nodes:
    node.SetSolutionStepValue(PATCH_INDEX, 1)
    node.Set(INTERFACE, True)
    
for i in range(0,int(len(elements)/4)):
    sys.stdout.flush()
    mesh1_model_part.CreateNewCondition("Condition3D4N", i+1,[int(elements[4*i+0]),int(elements[4*i+1]),int(elements[4*i+2]),int(elements[4*i+3])], mesh1_model_part.GetProperties()[1])

#for cond in mesh1_model_part.Conditions:
#    cond.Set(INTERFACE, True)
for cond in mesh1_model_part.Conditions:
    if cond.GetNode(0).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(1).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(2).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(3).GetSolutionStepValue(PATCH_INDEX) == 1.0 :
        cond.Set(INTERFACE, True)
    else:
        cond.Set(INTERFACE, False)


# Mesh 2 Definition
nodes = np.array([0.0,0.0,0.0, 1.0,0.0,0.0, 2.0,0.0,0.0, 3.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, 2.0,1.0,0.0, 3.0,1.0,0.0])
elements = np.array([1,2,6,5, 2,3,7,6, 3,4,8,7])

#nodes = np.array([0.0,0.0,0.0, 1.5,0.0,0.0, 3.0,0.0,0.0, 0.0,1.0,0.0, 1.5,1.0,0.0, 3.0,1.0,0.0])
#elements = np.array([1,2,5,4, 2,3,6,5])

#Adding model part for second Mesh
print("Mesh 2 ::::: Length of Nodes = ",int(len(nodes)/3))
print("Mesh 2 ::::: Length of Elements = ",int(len(elements)/4))
for i in range(0,int(len(nodes)/3)):
    mesh2_model_part.CreateNewNode(i+1, nodes[3*i+0], nodes[3*i+1], nodes[3*i+2])

for node in mesh2_model_part.Nodes:
    node.SetSolutionStepValue(PATCH_INDEX, 1)
    node.Set(INTERFACE, True)

mesh2_model_part.AddProperties(Properties(1))
for i in range(0,int(len(elements)/4)):
    sys.stdout.flush()
    mesh2_model_part.CreateNewCondition("Condition3D4N", i+1,[int(elements[4*i+0]),int(elements[4*i+1]),int(elements[4*i+2]),int(elements[4*i+3])], mesh2_model_part.GetProperties()[1])

for cond in mesh2_model_part.Conditions:
    if cond.GetNode(0).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(1).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(2).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(3).GetSolutionStepValue(PATCH_INDEX) == 1.0 :
        cond.Set(INTERFACE, True)
    else:
        cond.Set(INTERFACE, False)

## Creating field on mesh1 
## This is just a pseudo field 
## Field[0] = sin(x+y^2)
## Field[1] = cos(y+z^2)
## Field[2] = sin(z+x^2)
for node in mesh1_model_part.Nodes:
    x = node.X
    y = node.Y
    z = node.Z
    fx = np.sin(x + y*y)
    fy = np.cos(y + z*z)
    fz = np.sin(z + x*x)

    node.SetSolutionStepValue(DISPLACEMENT_X, fx)
    node.SetSolutionStepValue(DISPLACEMENT_Y, fy)
    node.SetSolutionStepValue(DISPLACEMENT_Z, fz)

mpr = RBFMapper.RBFMapper(mesh1_model_part,mesh2_model_part,2.0)
print(mpr.MapFromMasterToSlave(DISPLACEMENT, VELOCITY))

