import KratosMultiphysics
from KratosMultiphysics.DEMApplication.sp_2d_rigid_fem_algorithm import DEMAnalysisStage2DSpRigidFem

model = KratosMultiphysics.Model()
with open("ProjectParametersDEM.json",'r') as parameter_file:
    project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
DEMAnalysisStage2DSpRigidFem(model, project_parameters).Run()
