from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as KratosIGA
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTesselation(KratosUnittest.TestCase):
    def test_basic_tesselation(self):

        
        model = KratosMultiphysics.Model()
        iga_model_part =  model.CreateModelPart("IgaModelPart")

        

        with open("geometry.json",'r') as geometry_file:
            iga_geometry_parameters = KratosMultiphysics.Parameters( geometry_file.read())

        iga_geometry_reader = KratosIGA.BrepJsonIO()

        embedded_iga_modeler = KratosIGA.EmbeddedIgaModeler(iga_model_part)
        embedded_iga_modeler.ImportGeometry(iga_geometry_reader, iga_geometry_parameters)
        
        reference_results = [
            [7.96505190663491, 5.0], 
            [6.871642978763861, 3.832994984245571], 
            [6.22799531912592, 3.2046354306827642], 
            [5.681904794197063, 2.7428147747867073], 
            [5.072267373364987, 2.3055491352480537], 
            [4.558284195802144, 1.9952216530689135], 
            [4.002699142060996, 1.712110978432267], 
            [3.406153099390143, 1.4600332572954393], 
            [2.7713553255020362, 1.2413776954420173], 
            [2.103508001461179, 1.0568356701512598], 
            [1.4107307845723192, 0.9051298418674958], 
            [0.70448536126865, 0.7827432658695289], 
            [0.0, 0.6836485039400402], 
            [0.0,0.0], 
            [10.0, 0.0], 
            [10.0, 5.0]]

        current_result = embedded_iga_modeler.PrintParametricTessellation()

        self.assertEqual(len(current_result), len(reference_results))

        for reference_result, current_result in zip(reference_results, current_result):
            self.assertAlmostEqual(reference_result, current_result)


if __name__ == '__main__':
    KratosUnittest.main()
