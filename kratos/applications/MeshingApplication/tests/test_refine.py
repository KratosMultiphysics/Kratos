from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication
import math

def GetFilePath(fileName):
    import os
    return os.path.dirname(os.path.realpath(__file__)) + "/../../../kratos/tests/" + fileName

class TestRedistance(KratosUnittest.TestCase):

    def _ComputeVolume(self,Elements):
        vol = 0.0
        for elem in Elements:
            vol += elem.GetArea()
        return vol
    
    def _ComputeSurfaceArea(self,Conditions):
        area = 0.0
        for cond in Conditions:
            area += cond.GetArea()
        return area       
        
    def test_refine_all(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        KratosMultiphysics.ModelPartIO(GetFilePath("coarse_sphere")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)
        print(model_part)
        
        original_vol = self._ComputeVolume(model_part.Elements)
        #original_surface = self._ComputeSurfaceArea(model_part.Conditions)
        
        for elem in model_part.Elements:
            elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)
            
        refiner = KratosMultiphysics.MeshingApplication.LocalRefineTetrahedraMesh(model_part)
        
        refine_on_reference = False
        interpolate_internal_variables = True
        refiner.LocalRefineMesh( refine_on_reference, interpolate_internal_variables)
        
        refined_vol = self._ComputeVolume(model_part.Elements)
        refined_surface = self._ComputeSurfaceArea(model_part.Conditions)

        print(refined_surface)
        self.assertAlmostEqual(refined_vol, original_vol, 9)
        #self.assertAlmostEqual(refined_surface, original_surface, 9)
        
    def test_refine_half(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        KratosMultiphysics.ModelPartIO(GetFilePath("coarse_sphere")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)
        print(model_part)
        
        original_vol = self._ComputeVolume(model_part.Elements)
        #original_surface = self._ComputeSurfaceArea(model_part.Conditions)
        
        for elem in model_part.Elements:
            if(elem.Id % 2 == 0):
                elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)
            
        refiner = KratosMultiphysics.MeshingApplication.LocalRefineTetrahedraMesh(model_part)
        
        refine_on_reference = False
        interpolate_internal_variables = True
        refiner.LocalRefineMesh( refine_on_reference, interpolate_internal_variables)
        
        refined_vol = self._ComputeVolume(model_part.Elements)
        refined_surface = self._ComputeSurfaceArea(model_part.Conditions)

        print(refined_surface)
        self.assertAlmostEqual(refined_vol, original_vol, 9)
        #self.assertAlmostEqual(refined_surface, original_surface, 9)


if __name__ == '__main__':
    KratosUnittest.main()
