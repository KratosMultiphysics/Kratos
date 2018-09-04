from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_mechanics_analysis
import structural_mechanics_analysis_nonlinear_sensitivity


# begin first analysis ****************************************************************************************************

with open("PrimalParameters.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())

model_part_name_primal = ProjectParametersPrimal["problem_data"]["model_part_name"].GetString()
model_truss_primal = Model()

primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_truss_primal, ProjectParametersPrimal)

primal_analysis.Run()

# end first analysis ****************************************************************************************************

print("Primal analysis finished. Start now to compute curvatures.")

# begin second analysis ****************************************************************************************************

with open("PrimalParameters_curvature.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())

model_part_name_primal = ProjectParametersPrimal["problem_data"]["model_part_name"].GetString()

curvature_analysis = structural_mechanics_analysis_nonlinear_sensitivity.StructuralMechanicsAnalysisNLSensitivity(model_truss_primal, ProjectParametersPrimal)

curvature_analysis.Run()

# end second analysis ****************************************************************************************************

#model_part = model_truss_primal.GetModelPart(model_part_name_primal)
#print("")
#print("********************************************")
#print("")
#for node in model_part.Nodes:
#    print("curvature x = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_X))
#    print("curvature y = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_Y))
#    print("curvature z = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_Z))
#    print("")
#    print("first order sen x = ", node.GetValue(NL_SENSITIVITY_FIRST_ORDER_X))
#    print("first order sen y = ", node.GetValue(NL_SENSITIVITY_FIRST_ORDER_Y))
#    print("first order sen z = ", node.GetValue(NL_SENSITIVITY_FIRST_ORDER_Z))
#    print("")
#    print("second order sen x = ", node.GetValue(NL_SENSITIVITY_SECOND_ORDER_X))
#    print("second order sen y = ", node.GetValue(NL_SENSITIVITY_SECOND_ORDER_Y))
#    print("second order sen z = ", node.GetValue(NL_SENSITIVITY_SECOND_ORDER_Z))
#    print("")
#    print("********************************************")
#    print("")



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




