import KratosMultiphysics
import KratosMultiphysics.IgaApplication as Iga
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os
print(os.getpid())

class Iga5pElementTests(KratosUnittest.TestCase):
    def _check_results(self,node,displacement_results, rotation_results):
        # check that the results are exact on the node
        disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[1], displacement_results[1], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)

    # # thick simply supported beam under dead load, geometrically non-linear computation
    # def test_nonlinear_beam_thick(self):
    #     current_directory = os.getcwd()
    #     os.chdir(current_directory + "/applications/IgaApplication/tests/nonlinear_beam_thick_p2_nCP22")

    #     with open("ProjectParameters_5p.json",'r') as parameter_file:
    #         parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
    #     model = KratosMultiphysics.Model()
    #     simulation = StructuralMechanicsAnalysis(model,parameters)
    #     simulation.Run()
        
    #     # deflections
    #     structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
    #     nodes = structural_analysis_model_part.Nodes
    #     displacement_results = [-0.728556500276984 , 6.683872400521456e-14 , -2.3189469880503393]
    #     rotation_results     = [-0.0010054918263754388 , -3.5289832037502245e-15 , 0.0]
    #     self._check_results(nodes[11],displacement_results, rotation_results)
    #     # internal forces
    #     elements = structural_analysis_model_part.Elements
    #     stress_cauchy_top_11 = elements[1].Calculate(Iga.STRESS_CAUCHY_TOP_11, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(stress_cauchy_top_11, 2.0619352112032185, 10)
    #     stress_cauchy_bottom_11 = elements[1].Calculate(Iga.STRESS_CAUCHY_BOTTOM_11, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(stress_cauchy_bottom_11, 2.0619352112032194, 10)
    #     internal_force_11 = elements[1].Calculate(Iga.INTERNAL_FORCE_11, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(internal_force_11, 7.1982308852259385, 10)
    #     shear_force_1 = elements[1].Calculate(Iga.SHEAR_FORCE_1, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(shear_force_1, -8.697908300327537, 10)        
        
    #     os.chdir(current_directory)

    # thick simply supported beam under dead load, geometrically linear computed
    def test_linear_beam_thick(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/linear_beam_thick_p4_nCP5")

        with open("shell_5p_project_parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.0 , 0.0 , 0.16458333333332775]
        rotation_results     = [-1.8010813997401652e-16 , 6.2727600891321345e-15 , 0.0]
        self._check_results(nodes[3],displacement_results, rotation_results)
        
        os.chdir(current_directory)

    # Scoredlis-Lo-Roof under dead load, geometrically linear computed
    def test_linear_scordelis(self):
        current_directory = os.getcwd()
        os.chdir(current_directory + "/linear_scordelis_p4_nCP10")

        with open("shell_5p_project_parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysis(model,parameters)
        simulation.Run()
        
        structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
        nodes = structural_analysis_model_part.Nodes
        displacement_results = [0.1625015226073027 , 0.0003449714905228694 , -0.3087376917471691]
        rotation_results     = [1.2016573629904487e-05 , 6.0947919423266035e-06 , 0.0]
        
        disp = nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        # The displacement in y direction is not checked because the model is kinematic in this direction and, therefore, the results can vary
        self.assertAlmostEqual(disp[0], displacement_results[0], 10)
        self.assertAlmostEqual(disp[2], displacement_results[2], 10)

        rot = nodes[5].GetSolutionStepValue(KratosMultiphysics.ROTATION)
        self.assertAlmostEqual(rot[0], rotation_results[0], 10)
        self.assertAlmostEqual(rot[1], rotation_results[1], 10)
        self.assertAlmostEqual(rot[2], rotation_results[2], 10)

        os.chdir(current_directory)

    # # thick simply supported beam under dead load, geometrically non-linear computed, review of computation of 
    # #   internal forces        
    # def test_nonlinear_beam_internal_forces(self):
    #     current_directory = os.getcwd()
    #     os.chdir(current_directory + "/applications/IgaApplication/tests/nonlinear_beam_thick_sd_p3_nCP83")

    #     with open("ProjectParameters_5p.json",'r') as parameter_file:
    #         parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
    #     model = KratosMultiphysics.Model()
    #     simulation = StructuralMechanicsAnalysis(model,parameters)
    #     simulation.Run()
        
    #     structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")
    #     elements = structural_analysis_model_part.Elements
    #     moment_result = elements[633].Calculate(Iga.INTERNAL_MOMENT_11, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(moment_result, -24.997442433906905, 10)
    #     shear_result = elements[1].Calculate(Iga.SHEAR_FORCE_1, structural_analysis_model_part.ProcessInfo)
    #     self.assertAlmostEqual(shear_result, -9.984730264779031, 10)

    #     os.chdir(current_directory)

if __name__ == '__main__':
    KratosUnittest.main()
