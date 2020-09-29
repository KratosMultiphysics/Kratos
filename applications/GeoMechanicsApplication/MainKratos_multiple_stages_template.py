from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
# import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication

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
