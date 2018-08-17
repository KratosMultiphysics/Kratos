from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_mechanics_analysis
import numpy as np
import matplotlib.pyplot as plt

def compute_EF_disp_curvature(model_part):
    load_factor = np.array([0.8, 0.9, 1.0])
    for node in model_part.Nodes:
        disp_1_x = node.GetSolutionStepValue(DISPLACEMENT_X, 3)
        disp_2_x = node.GetSolutionStepValue(DISPLACEMENT_X, 2)
        disp_3_x = node.GetSolutionStepValue(DISPLACEMENT_X, 1)
        disp_1_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 3)
        disp_2_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 2)
        disp_3_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 1)
        disp_1_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 3)
        disp_2_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 2)
        disp_3_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 1)

        disp_x = np.array([disp_1_x, disp_2_x, disp_3_x])
        px = np.polyfit(load_factor, disp_x, 2)
        disp_y = np.array([disp_1_y, disp_2_y, disp_3_y])
        py = np.polyfit(load_factor, disp_y, 2)
        disp_z = np.array([disp_1_z, disp_2_z, disp_3_z])
        pz = np.polyfit(load_factor, disp_z, 2)

        node.SetValue(DISPLACEMENT_NL_SENSITIVITY_X, 2 * px[0])
        node.SetValue(DISPLACEMENT_NL_SENSITIVITY_Y, 2 * py[0])
        node.SetValue(DISPLACEMENT_NL_SENSITIVITY_Z, 2 * pz[0])
        #print("curvature x = ", 2 * px[0])
        #print("curvature y = ", 2 * py[0])
        #print("curvature z = ", 2 * pz[0])
        #print("")


with open("PrimalParameters.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())

model_part_name_primal = ProjectParametersPrimal["problem_data"]["model_part_name"].GetString()
model_truss_primal = Model()

primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_truss_primal, ProjectParametersPrimal)

primal_analysis.Initialize()

model_part = model_truss_primal.GetModelPart(model_part_name_primal)
model_part.SetBufferSize(5)

primal_analysis.RunSolutionLoop()
compute_EF_disp_curvature(model_part)

primal_analysis.Finalize()

print("")

for node in model_part.Nodes:
    print("curvature x = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_X))
    print("curvature y = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_Y))
    print("curvature z = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_Z))
    print("")



#print("")
#print("polynomial coeff = ", z)
#print("curvature = ", 2 * z[0])
#print("")
#print(disp_1)
#print(disp_2)
#print(disp_3)
#print("Finished polynomial interpolation!")

#p = np.poly1d(z)
#print(np.poly1d(p))
#xp = np.linspace(0, 1.0, 100)
#_ = plt.plot(x, y, '.', xp, p(xp), '-')
#plt.ylim(0,0.2)
#plt.show()




