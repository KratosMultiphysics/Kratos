from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.IgaApplication as Iga
from KratosMultiphysics.IgaApplication.iga_analysis import IgaAnalysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os
print(os.getpid())

class Iga7pElementTests(KratosUnittest.TestCase):
    def _check_results(self,node,displacement_results, rotation_results, w_bar_result):
        #check that the results are exact on the node
        disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[1], displacement_results[1], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)

        w_bar = node.GetSolutionStepValue(Iga.W_BAR)
        self.assertAlmostEqual(w_bar, w_bar_result, 10)

    def test_linear_beam_thick(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/linear_beam_thick_p4_nCP5")

        with open("ProjectParameters_7p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.0 , 0.0 , 0.16458333333332836]
        rotation_results     = [-8.554774384372251e-17, -1.0923034942261134e-18 , 0.0]
        w_bar_result = 0.0
        self._check_results(nodes[3],displacement_results, rotation_results, w_bar_result)
        
        os.chdir(current_directory)

    def test_linear_plate(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/linear_plate_Lt5_p2_nCP12")

        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.0 , 0.0 , -0.5837065902297612]
        rotation_results     = [-0.002635355333859876 , -0.002635355333773054 , 0.0]
        w_bar_result = 0.019122079777579713
        self._check_results(nodes[66],displacement_results, rotation_results, w_bar_result)
        
        os.chdir(current_directory)

    def test_linear_scordelis(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/applications/IgaApplication/tests/linear_scordelis_p4_nCP10")

        with open("ProjectParameters_7p.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = IgaAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.16250201185879148 , 0.06159704287645379 , -0.30873854692330444]
        rotation_results     = [4.6836363064969345e-06 , 1.2016618499442812e-05 , 3.808622928110629e-06]
        w_bar_result = -8.37565509156725e-10

        disp = nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        # The displacement in y direction is not checked because the model is kinematic in this direction and, therefore, the results can vary
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = nodes[5].GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)

        w_bar = nodes[5].GetSolutionStepValue(Iga.W_BAR)
        self.assertAlmostEqual(w_bar, w_bar_result, 10)

        os.chdir(current_directory)
            
    # def test_nonlinear_beam_internal_forces(self):
    #     current_directory = os.getcwd()
    #     os.chdir(current_directory + "/applications/IgaApplication/tests/nonlinear_beam_thick_sd_p3_nCP83")

    #     with open("ProjectParameters_5p.json",'r') as parameter_file:
    #         parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
    #     model = KratosMultiphysics.Model()
    #     simulation = IgaAnalysis(model,parameters)
    #     simulation.Run()
        
    #     structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
    #     elements = structural_analysis_model_part.Elements
    #     moment_result = elements[633].Calculate(Iga.INTERNAL_MOMENT_11, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(moment_result, -24.997443146162688, 10)
    #     shear_result = elements[1].Calculate(Iga.SHEAR_FORCE_1, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(shear_result, -9.984730264779031, 10)

    #     os.chdir(current_directory)

if __name__ == '__main__':
    KratosUnittest.main()