import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.FemToDemApplication
import CouplingFemDem3D
import CouplingFemDem3DHexahedrons
	
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())
problem_data = ProjectParameters["problem_data"]
model = KratosMultiphysics.Model()

if problem_data.Has("is_hexahedron") and problem_data["is_hexahedron"].GetBool() == True:
	CouplingFemDem3DHexahedrons.FEMDEM3DHexahedrons_Solution(model).Run()	
else:
	CouplingFemDem3D.FEMDEM3D_Solution(model).Run()