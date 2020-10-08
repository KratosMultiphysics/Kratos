import os
import KratosMultiphysics

from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

if __name__ == "__main__":

	project_paths = [r"path\first_stage",
					r"path\second_stage",
					r"path\third_stage",
					r"path\fourth_stage",
	]

	model = KratosMultiphysics.Model()
	for project_path in project_paths: 
		os.chdir(project_path)
		
		with open("ProjectParameters.json",'r') as parameter_file:
			parameters = KratosMultiphysics.Parameters(parameter_file.read())

		simulation = GeoMechanicsAnalysis(model,parameters)
		simulation.Run()
