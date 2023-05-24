import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_analysis import HelmholtzAnalysis
from KratosMultiphysics.testing.utilities import ReadModelPart
import KratosMultiphysics.OptimizationApplication as KOA
from numpy import linalg as LA
import os


model = KM.Model()
model_part = model.CreateModelPart("hl_solid_shell")

model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VECTOR)
KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_X, model_part)
KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Y, model_part)
KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VECTOR_Z, model_part)

model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VARS_SHAPE)
model_part.AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SOURCE_SHAPE)
KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VARS_SHAPE_X, model_part)
KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VARS_SHAPE_Y, model_part)
KM.VariableUtils().AddDof(KOA.HELMHOLTZ_VARS_SHAPE_Z, model_part)

model_part.CreateNewNode(1, 0.0,0.0,0.0)
model_part.CreateNewNode(2, 1.0,0.0,0.0)
model_part.CreateNewNode(3, 0.0,1.0,0.0)
model_part.CreateNewNode(4, 0.0,0.0,1.0)

model_part.CreateNewNode(5, 2.0,0.0,0.0)
model_part.CreateNewNode(6, 3.0,0.0,0.0)
model_part.CreateNewNode(7, 2.0,1.0,0.0)
model_part.CreateNewNode(8, 2.0,0.0,1.0)

for node in model_part.Nodes:
    node.SetValue(KM.NUMBER_OF_NEIGHBOUR_ELEMENTS,1)
    node.SetValue(KOA.HELMHOLTZ_VECTOR_SOURCE,[1,1,1])
    node.SetSolutionStepValue(KOA.HELMHOLTZ_SOURCE_SHAPE,[1,1,1])


properties_1 = model_part.CreateNewProperties(1)
properties_2 = model_part.CreateNewProperties(2)
properties_2.SetValue(KOA.HELMHOLTZ_BULK_RADIUS_SHAPE, 1000.0)
properties_2.SetValue(KOA.HELMHOLTZ_SURF_RADIUS_SHAPE, 5.0)
model_part.ProcessInfo.SetValue(KOA.COMPUTE_HELMHOLTZ_INVERSE, False)
model_part.ProcessInfo.SetValue(KOA.COMPUTE_CONTROL_POINTS_SHAPE, False)
model_part.ProcessInfo.SetValue(KOA.HELMHOLTZ_INTEGRATED_FIELD, True)
model_part.ProcessInfo.SetValue(KOA.HELMHOLTZ_BULK_RADIUS_SHAPE, 1000.0)
model_part.ProcessInfo.SetValue(KOA.HELMHOLTZ_RADIUS, 5.0)


model_part.CreateNewElement("HelmholtzSolidShapeElement3D4N", 1, [1,2,3,4], properties_1)
model_part.CreateNewCondition("HelmholtzSurfaceShapeCondition3D3N", 1, [2,1,3], properties_1)
model_part.CreateNewElement("HelmholtzBulkShape3D4N", 2, [5,6,7,8], properties_2)
model_part.CreateNewCondition("HelmholtzSurfShapeCondition3D3N", 2, [6,5,7], properties_2)

tmoc = KM.TetrahedralMeshOrientationCheck
flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
KM.TetrahedralMeshOrientationCheck(model_part, False, flags).Execute()

LHS_1 = KM.Matrix()
RHS_1 = KM.Vector()

#model_part.GetElement(1).CalculateLocalSystem(LHS_1, RHS_1, model_part.ProcessInfo)
##print("LHS_1: ",LHS_1)
#test_strain_1 = model_part.GetElement(1).Calculate(KOA.ELEMENT_STRAIN_ENERGY,model_part.ProcessInfo)
#print("test_strain_1: ",test_strain_1)
##print("RHS_1: ",RHS_1)


#LHS_2 = KM.Matrix()
#RHS_2 = KM.Vector()
#model_part.GetElement(2).CalculateLocalSystem(LHS_2, RHS_2, model_part.ProcessInfo)
#test_strain_2 = model_part.GetElement(2).Calculate(KOA.ELEMENT_STRAIN_ENERGY,model_part.ProcessInfo)
#print("test_strain_2: ",test_strain_2)
##print("LHS_1-LHS_2: ",LHS_1-LHS_2)
##print("RHS_2: ",RHS_2)

model_part.ProcessInfo.SetValue(KOA.COMPUTE_HELMHOLTZ_INVERSE, False)

c_LHS_1 = KM.Matrix()
c_RHS_1 = KM.Vector()
model_part.GetCondition(1).CalculateLocalSystem(c_LHS_1, c_RHS_1, model_part.ProcessInfo)
c_test_strain_1 = model_part.GetCondition(1).Calculate(KOA.ELEMENT_STRAIN_ENERGY,model_part.ProcessInfo)
print("c_test_strain_1: ",c_test_strain_1)
print("c_LHS_1: ",c_LHS_1)
print("c_RHS_1: ",c_RHS_1)

#c_LHS_2 = KM.Matrix()
#c_RHS_2 = KM.Vector()
#model_part.GetCondition(2).CalculateLocalSystem(c_LHS_2, c_RHS_2, model_part.ProcessInfo)
#c_test_strain_2 = model_part.GetCondition(2).Calculate(KOA.ELEMENT_STRAIN_ENERGY,model_part.ProcessInfo)
#print("c_test_strain_2: ",c_test_strain_2)
#print("c_LHS_2: ",c_LHS_2)
#print("c_RHS_2: ",c_RHS_2)
#print("c_LHS_1-c_LHS_2: ",c_LHS_1-c_LHS_2)
