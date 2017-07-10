from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
import numpy as np
from numpy import linalg as la
#
#
# setting the domain size for the problem to be solved
#
import sys
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.EmpireApplication import *
import empire_mapper

# defining variables to be used
# GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
# CODE IS NOT USING IT

variables_dictionary = {"DISPLACEMENT": DISPLACEMENT,
                        "PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "REACTION": REACTION,
                        "DISTANCE": DISTANCE, }

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
print("Mesh 1 ::::: Maximum entry in elements vector = ", np.max(elements))
for i in range(0,int(len(nodes)/3)):
    mesh1_model_part.CreateNewNode(i+1, nodes[3*i+0], nodes[3*i+1], nodes[3*i+2])
print("Mesh 1 ::::: Finished adding nodes to model part !!")

# Creating elements with the above nodes
mesh1_model_part.AddProperties(Properties(1))
print("Mesh 1 ::::: Finished adding Properties to model part !!")
prp = mesh1_model_part.GetProperties()[1]
print("Mesh 1 ::::: Finished extracting properties from the model part !!")

for node in mesh1_model_part.Nodes:
    node.SetSolutionStepValue(PATCH_INDEX, 1)
    node.Set(INTERFACE, True)
    
for i in range(0,int(len(elements)/4)):
    sys.stdout.flush()
    mesh1_model_part.CreateNewCondition("Condition3D4N", i+1,[int(elements[4*i+0]),int(elements[4*i+1]),int(elements[4*i+2]),int(elements[4*i+3])], mesh1_model_part.GetProperties()[1])
print("Mesh 1 ::::: Finished Adding conditions to model part !!")

#for cond in mesh1_model_part.Conditions:
#    cond.Set(INTERFACE, True)
print("Mesh 1 ::::: Finished setting conditions in model part !!")
for cond in mesh1_model_part.Conditions:
    if cond.GetNode(0).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(1).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(2).GetSolutionStepValue(PATCH_INDEX) == 1.0 \
       and cond.GetNode(3).GetSolutionStepValue(PATCH_INDEX) == 1.0 :
        cond.Set(INTERFACE, True)
    else:
        cond.Set(INTERFACE, False)

print("Mesh 1 ::::: Finished adding Properties to model part !! Number of the conditions added to model part are : ", len(mesh1_model_part.Conditions))



# Mesh 2 Definition
nodes = np.array([0.0,0.0,0.0, 1.0,0.0,0.0, 2.0,0.0,0.0, 3.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, 2.0,1.0,0.0, 3.0,1.0,0.0])
elements = np.array([1,2,6,5, 2,3,7,6, 3,4,8,7])

#nodes = np.array([0.0,0.0,0.0, 1.5,0.0,0.0, 3.0,0.0,0.0, 0.0,1.0,0.0, 1.5,1.0,0.0, 3.0,1.0,0.0])
#elements = np.array([1,2,5,4, 2,3,6,5])

#Adding model part for second Mesh
print("Mesh 2 ::::: Length of Nodes = ",int(len(nodes)/3))
print("Mesh 2 ::::: Length of Elements = ",int(len(elements)/4))
print("Mesh 2 ::::: Maximum entry in elements vector = ", np.max(elements))
for i in range(0,int(len(nodes)/3)):
    mesh2_model_part.CreateNewNode(i+1, nodes[3*i+0], nodes[3*i+1], nodes[3*i+2])

for node in mesh2_model_part.Nodes:
    node.SetSolutionStepValue(PATCH_INDEX, 1)
    node.Set(INTERFACE, True)

mesh2_model_part.AddProperties(Properties(1))
for i in range(0,int(len(elements)/4)):
    sys.stdout.flush()
    mesh2_model_part.CreateNewCondition("Condition3D4N", i+1,[int(elements[4*i+0]),int(elements[4*i+1]),int(elements[4*i+2]),int(elements[4*i+3])], mesh2_model_part.GetProperties()[1])

print("Mesh 2 ::::: Finished setting conditions in model part !!")
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


## Creating the mapper
isDual = 0
isOppSurfNormal = 0
isEnforceConsistancy = 1
mapper = empire_mapper.FiniteElementMapper("default",3,mesh1_model_part,mesh2_model_part,isDual,isOppSurfNormal,isEnforceConsistancy)
print("############### Successfully created mapper !!")

## Mapping the DISPLACEMENT field of mesh1_model_part  to VELOCITY field of mesh2_model_part
## The fields DISPLACEMENT and VELOCITY here are pseudo and not physical.
mapper.doMapping(mesh1_model_part,DISPLACEMENT,mesh2_model_part,VELOCITY)

norm = 0
i = 0
## Validating the mapped values
for node in mesh2_model_part.Nodes:
    cvalx = node.GetSolutionStepValue(VELOCITY_X) # cval* is calculated value 
    cvaly = node.GetSolutionStepValue(VELOCITY_Y)
    cvalz = node.GetSolutionStepValue(VELOCITY_Z)
    #
    # Calculating the analytical value
    #
    x = node.X
    y = node.Y
    z = node.Z
    avalx = np.sin(x + y*y) # aval* is the analytical value 
    avaly = np.cos(y + z*z)
    avalz = np.sin(z + x*x)
    print('Node :: :: :: ', i)
    print('Analytical x :: ',avalx, ' :: analytical y :: ', avaly, ' :: analytical z :: ',avalz)
    print('calculated x :: ',cvalx, ' :: calculated y :: ', cvaly, ' :: calculated z :: ',cvalz)
    
    norm += la.norm(np.array([cvalx-avalx, cvaly-avaly, cvalz-avalz]), np.inf)
    i += 1
    
print("Norm of the difference between the mapped and analytical fields is : ",norm)

    
