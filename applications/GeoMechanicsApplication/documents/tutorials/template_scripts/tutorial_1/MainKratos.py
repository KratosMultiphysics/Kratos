from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import KratosMultiphysics
import KratosMultiphysics.LinearSolversApplication
# import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication

from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

if __name__ == "__main__":

	project_paths = [r"D:\software_development\Kratos\applications\GeoMechanicsApplication\test_examples\tutorial_1_stage_1.gid",
					r"D:\software_development\Kratos\applications\GeoMechanicsApplication\test_examples\tutorial_1_stage_2.gid",
					r"D:\software_development\Kratos\applications\GeoMechanicsApplication\test_examples\tutorial_1_stage_3.gid",
					r"D:\software_development\Kratos\applications\GeoMechanicsApplication\test_examples\tutorial_1_stage_4.gid",
	]

	model = KratosMultiphysics.Model()
	for project_path in project_paths: 
		os.chdir(project_path)
		
		with open("ProjectParameters.json",'r') as parameter_file:
			parameters = KratosMultiphysics.Parameters(parameter_file.read())

		simulation = GeoMechanicsAnalysis(model,parameters)
		simulation.Run()
