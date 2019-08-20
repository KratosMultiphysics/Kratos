from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.IgaApplication.iga_analysis import IgaAnalysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os
print(os.getpid())

class Iga5pElementTests(KratosUnittest.TestCase):
    def _check_results(self,node,displacement_results, rotation_results):
        #check that the results are exact on the node
        disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[1], displacement_results[1], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)

    def test_nonlinear_beam_thick(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/nonlinear_beam_thick_p2_nCP22")

        with open("ProjectParameters_5p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [-0.728556500276984 , 6.683872400521456e-14 , -2.3189469880503393]
        rotation_results     = [-0.0010054918263754388 , -3.5289832037502245e-15 , 0.0]
        self._check_results(nodes[11],displacement_results, rotation_results)
        
        os.chdir(current_directory)

    def test_linear_beam_thick(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/linear_beam_thick_p4_nCP5")

        with open("ProjectParameters_5p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.0 , 0.0 , 0.16458333333332775]
        rotation_results     = [-1.8010813997401652e-16 , 6.2727600891321345e-15 , 0.0]
        self._check_results(nodes[3],displacement_results, rotation_results)
        
        os.chdir(current_directory)

    def test_linear_scordelis(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/linear_scordelis_p4_nCP10")

        with open("ProjectParameters_5p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.1625099984030851 , -0.5041700440919504 , -0.30870341394117473]
        rotation_results     = [1.202168535482754e-05 , 6.0867104344852356e-06 , 0.0]
        
        disp = nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        # The displacement in y direction is not checked because the model is kinematic in this direction and, therefore, the results can vary
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = nodes[5].GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)

        os.chdir(current_directory)

if __name__ == '__main__':
    KratosUnittest.main()