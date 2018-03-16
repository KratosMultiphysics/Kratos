from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import math, time
from KratosMultiphysics import *
#from KratosMultiphysics.AltairThermalApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
#from KratosMultiphysics.Click2CastApplication import *
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.SwimmingDEMApplication as dem


import KratosMultiphysics
#import KratosMultiphysics.AltairThermalApplication as KratosAltairThermal

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

#create model part
main_model_part = ModelPart("ModelPart")
main_model_part.AddProperties(Properties(1))

#add variables
main_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
main_model_part.AddNodalSolutionStepVariable(TEMPERATURE_GRADIENT)
main_model_part.AddNodalSolutionStepVariable(NODAL_VOLUME)
main_model_part.AddNodalSolutionStepVariable(VELOCITY)
main_model_part.AddNodalSolutionStepVariable(VELOCITY_Z_GRADIENT)



#populate model_part with nodes and elements
new_node_id = 1
new_elem_id = 1
n = 10
h = 1/n

for k in range(1,n+2):
   for j in range(1,n+2):
       for i in range(1,n+2):
           main_model_part.CreateNewNode(new_node_id,h*(i-1),h*(j-1),h*(k-1))
           new_node_id=new_node_id+1

# We do a loop on hexahedras and then we obtain 6 tetrahedon from each hexahedra

for k in range(0,n):
   for j in range(0,n):
       for i in range(0,n):

           #Find the 8 nodes that are present in the hexahedra
           node0 = (n+1)*(n+1)*k+(n+1)*j+i+1
           node1 = (n+1)*(n+1)*(k+1)+(n+1)*j+i+1
           node2 = (n+1)*(n+1)*(k+1)+(n+1)*j+i+2
           node3 = (n+1)*(n+1)*k+(n+1)*j+i+2
           node4 = (n+1)*(n+1)*k+(n+1)*(j+1)+i+1
           node5 = (n+1)*(n+1)*(k+1)+(n+1)*(j+1)+i+1
           node6 = (n+1)*(n+1)*(k+1)+(n+1)*(j+1)+i+2
           node7 = (n+1)*(n+1)*k+(n+1)*(j+1)+i+2

           #Using the 6 different combinations we obtain the 6 tetrahedon
           main_model_part.CreateNewElement("Element3D4N", new_elem_id, [node3,node2,node6,node0], main_model_part.GetProperties()[1])
           new_elem_id = new_elem_id+1
           main_model_part.CreateNewElement("Element3D4N", new_elem_id, [node0,node3,node7,node6], main_model_part.GetProperties()[1])
           new_elem_id = new_elem_id+1
           main_model_part.CreateNewElement("Element3D4N", new_elem_id, [node1,node5,node6,node0], main_model_part.GetProperties()[1])
           new_elem_id = new_elem_id+1
           main_model_part.CreateNewElement("Element3D4N", new_elem_id, [node0,node4,node5,node6], main_model_part.GetProperties()[1])
           new_elem_id = new_elem_id+1
           main_model_part.CreateNewElement("Element3D4N", new_elem_id, [node0,node1,node2,node6], main_model_part.GetProperties()[1])
           new_elem_id = new_elem_id+1
           main_model_part.CreateNewElement("Element3D4N", new_elem_id, [node4,node7,node6,node0], main_model_part.GetProperties()[1])
           new_elem_id = new_elem_id+1


print("number of nodes", main_model_part.NumberOfNodes())
print("number of elements", main_model_part.NumberOfElements())

for node in main_model_part.Nodes:
    node.SetValue(core.NODAL_VOLUME, 0.0)
    for element in main_model_part.Elements:
        el_vol = element.GetArea()
        for node in element.GetNodes():
            nodal_volume = node.GetValue(core.NODAL_VOLUME) + el_vol*0.25
            node.SetValue(core.NODAL_VOLUME, nodal_volume)

#add dofs
for node in main_model_part.Nodes:
    node.AddDof(KratosMultiphysics.TEMPERATURE_GRADIENT_X)
    node.AddDof(KratosMultiphysics.TEMPERATURE_GRADIENT_Y)
    node.AddDof(KratosMultiphysics.TEMPERATURE_GRADIENT_Z)



K_diff = 0.00001         #thermal diffusivity = (conductivity/(density*specific_heat))
T0 = 1000          #initial temperature
t = 1               #time
pi = math.pi

# model_part_io = ModelPartIO(GetFilePath("Levers_2"))
# model_part_io.ReadModelPart(main_model_part)

#we introduce a Fourier series expansion for a the Temperature Field T(x,y,z) = 4 * T0/pi * sin(pi x) * exp(-k*pi^2*t)
for node in main_model_part.Nodes:
    u = 0
    #m = 127           # number of terms of fourier seres = (2*m-1)
    #for i in range (1,m+1):
    #    k = 2*i - 1
    #    b = 4*T0/((2*k-1)*pi)
    #   u += b*math.sin(k*pi*node.X)*math.exp(-K_diff*pi*pi*k*k*t)
    u = math.sin(pi*node.X)*math.sin(pi*node.Y)*math.sin(pi*node.Z)

    # u = -2*node.X + 3*node.Y + 3*node.Z
    # u = math.exp(-25*node.X) + math.exp(-25*node.Y) + math.exp(-25*node.Z)
    #u = math.sinh(node.X)
    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, u)
    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_X, 0)
    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_Y, 0)
    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_Z, 0)

#Construct and Execute the process that calculates the temperature gradient
niter_max = 1

#compute_gradient1 = KratosAltairThermal.ReconstructNodalGradientProcess(main_model_part , KratosMultiphysics.TEMPERATURE , KratosMultiphysics.TEMPERATURE_GRADIENT, niter_max)
#compute_gradient.Execute()
import derivative_recover_process
compute_gradient=derivative_recover_process.CreateProcess(main_model_part, KratosMultiphysics.TEMPERATURE , KratosMultiphysics.TEMPERATURE_GRADIENT)

start = time.time()
compute_gradient.Execute()
end = time.time()
final_time = str(end - start)


#Absolute and relative errors
abs_err = 0
rel_err = 0
max_abs_err = 0
max_rel_err = 0
error_l2 = 0
for node in main_model_part.Nodes:
    Tx = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_X)
    Ty = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_Y)
    Tz = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_Z)

    # Tx_ana = -25*math.exp(-25*node.X)
    # Ty_ana = -25*math.exp(-25*node.Y)
    # Tz_ana = -25*math.exp(-25*node.Z)
    Tx_ana = math.cos(pi*node.X)*math.sin(pi*node.Y)*math.sin(pi*node.Z)*pi
    Ty_ana = math.sin(pi*node.X)*math.cos(pi*node.Y)*math.sin(pi*node.Z)*pi
    Tz_ana = math.sin(pi*node.X)*math.sin(pi*node.Y)*math.cos(pi*node.Z)*pi

    #error_l2 += node.GetValue(KratosMultiphysics.NODAL_VOLUME)*((Tx-Tx_ana)**2 + (Ty-Ty_ana)**2 + (Tz-Tz_ana)**2)
    error_l2 += node.GetValue(NODAL_VOLUME)*((Tx-Tx_ana)**2 + (Ty-Ty_ana)**2 + (Tz-Tz_ana)**2)


error_l2 = error_l2**(1/2)

    #Analitical derivative of the temperature field
#    Tx_ana = 0
#    for i in range (1,m+1):
#        k = 2*i - 1
#        b = 4*T0/(2*k-1)
#        Tx_ana += b*math.cos(k*pi*node.X)*math.exp(-K_diff*pi*pi*k*k*t)*k


#    abs_err_node = abs(Tx - Tx_ana) + abs(Ty) + abs(Tz)

#    abs_err += abs_err_node
#    err_rel_node = abs_err_node/Tx_ana
#    rel_err += err_rel_node

#    #print(Tx_ana)
#    #print(node.GetSolutionStepValue(TEMPERATURE_GRADIENT_X))
#    #print(err_rel_node)

#    if abs_err_node > max_abs_err:
#        max_abs_err = abs_err_node
#    if err_rel_node > max_rel_err:
#        max_rel_err = err_rel_node
#abs_err /= main_model_part.NumberOfNodes()
#rel_err /= main_model_part.NumberOfNodes()

#print("Average abs error", abs_err)
#print("Maximum abs error", max_abs_err)
#print("Average rel error", rel_err)
#print("Maximum rel error", max_rel_err)

#write results

target = open("solution.txt", 'w')
target.truncate()
# for node in main_model_part.Nodes:
#     target.write(str(node.Id)+" " + str(node.X0)+" " + str(node.Y0)+" "+ str(node.Z0)+" " +str(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z))+" " +str(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z_GRADIENT_X))+" " +str(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z_GRADIENT_Y))+" " +str(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z_GRADIENT_Z)))
# #     target.write(str(node.X) + "   " + str(-25*math.exp(-25*node.X)) + "   " + str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_X))+ "   " + str(abs(math.cosh(node.X) - node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT_X))) )
#     target.write("\n")

target.write(str(error_l2))
target.write("\n")
target.write(str(final_time))

target.close()


