from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.IgaApplication as Iga
from KratosMultiphysics.IgaApplication.iga_analysis import IgaAnalysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os
print(os.getpid())

class Iga3pElementTests(KratosUnittest.TestCase):
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

        with open("ProjectParameters_3p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [-0.7140554372203562 , -3.795961212461202e-16 , -2.298216150730901]
        rotation_results     = [0.0 , 0.0 , 0.0]
        self._check_results(nodes[11],displacement_results, rotation_results)

        os.chdir(current_directory)

    def test_linear_beam_thick(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/linear_beam_thick_p4_nCP5")

        with open("ProjectParameters_3p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.0 , 0.0 , 0.15624999999999398]
        rotation_results     = [0.0 , 0.0 , 0.0]
        self._check_results(nodes[3],displacement_results, rotation_results)
        
        os.chdir(current_directory)

    def test_linear_scordelis(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/linear_scordelis_p4_nCP10")

        with open("ProjectParameters_3p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.16237062769605737 , 0.05098425540534873 , -0.3084820432701192]
        rotation_results     = [0.0 , 0.0 , 0.0]

        disp = nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        # The displacement in y direction is not checked because the model is kinematic in this direction and, therefore, the results can vary
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = nodes[5].GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)

        os.chdir(current_directory)

    def test_nonlinear_beam_internal_forces(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/nonlinear_beam_thick_sd_p3_nCP83")

        with open("ProjectParameters_3p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        elements = structural_analysis_model_part.Elements
        moment_result = elements[633].Calculate(Iga.INTERNAL_MOMENT_11, structural_analysis_model_part.ProcessInfo)
        self.assertAlmostEqual(moment_result, -24.99744280630455, 10)
        shear_result = elements[1].Calculate(Iga.SHEAR_FORCE_1, structural_analysis_model_part.ProcessInfo)
        self.assertAlmostEqual(shear_result, -9.874997238324553, 10)

        os.chdir(current_directory)

if __name__ == '__main__':
    KratosUnittest.main()