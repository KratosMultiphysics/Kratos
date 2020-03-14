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
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        KratosMultiphysics.ModelPartIO(GetFilePath("coarse_sphere")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        #IMPORTANT! must compute the neighbours first!!
        KratosMultiphysics.FindNodalNeighboursProcess(model_part).Execute()
        
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

        print("refined_vol = ",refined_vol)
        self.assertAlmostEqual(refined_vol, original_vol, 12)
        #self.assertAlmostEqual(refined_surface, original_surface, 9)
        self.assertEqual(len(model_part.Elements), 1992)
        self.assertEqual(len(model_part.Nodes), 482)
        
    def test_refine_half(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        KratosMultiphysics.ModelPartIO(GetFilePath("coarse_sphere")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)
        print(model_part)
        
        original_vol = self._ComputeVolume(model_part.Elements)
        #original_surface = self._ComputeSurfaceArea(model_part.Conditions)
        
        for elem in model_part.Elements:
            if(elem.Id % 2 == 0):
                elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)
            
        #IMPORTANT! must compute the neighbours first!!
        KratosMultiphysics.FindNodalNeighboursProcess(model_part).Execute()
        
        refiner = KratosMultiphysics.MeshingApplication.LocalRefineTetrahedraMesh(model_part)
        
        refine_on_reference = False
        interpolate_internal_variables = True
        refiner.LocalRefineMesh( refine_on_reference, interpolate_internal_variables)
        
        refined_vol = self._ComputeVolume(model_part.Elements)
        refined_surface = self._ComputeSurfaceArea(model_part.Conditions)

        print("refined_vol = ",refined_vol)
        self.assertAlmostEqual(refined_vol, original_vol, 12)
        #self.assertAlmostEqual(refined_surface, original_surface, 9)
        self.assertEqual(len(model_part.Elements), 2010)
        self.assertEqual(len(model_part.Nodes), 462)
                         
    def test_refine_half_and_improve(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        KratosMultiphysics.ModelPartIO(GetFilePath("coarse_sphere")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)
        print(model_part)
        
        original_vol = self._ComputeVolume(model_part.Elements)
        #original_surface = self._ComputeSurfaceArea(model_part.Conditions)
        
        for elem in model_part.Elements:
            if(elem.Id % 2 == 0):
                elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)
            
        #IMPORTANT! must compute the neighbours first!!
        KratosMultiphysics.FindNodalNeighboursProcess(model_part).Execute()
        
        refiner = KratosMultiphysics.MeshingApplication.LocalRefineTetrahedraMesh(model_part)
        
        refine_on_reference = False
        interpolate_internal_variables = True
        refiner.LocalRefineMesh( refine_on_reference, interpolate_internal_variables)

        KratosMultiphysics.FindNodalNeighboursProcess(model_part).Execute()
        
        reconnector = KratosMultiphysics.MeshingApplication.TetrahedraReconnectUtility(model_part)
        simIter = 2
        iterations = 2
        ProcessByNode = False
        ProcessByFace = True
        ProcessByEdge = True
        saveToFile = False
        removeFreeVertexes = False
        evaluateInParallel = True
        reInsertNodes = False
        debugMode = False
        minAngle = 15


        reconnector.setBlockSize(2048)
        reconnector.OptimizeQuality(model_part, simIter, iterations, ProcessByNode, ProcessByFace, ProcessByEdge, saveToFile, removeFreeVertexes, evaluateInParallel, reInsertNodes, debugMode,minAngle)
        meshIsValid = reconnector.EvaluateQuality()
        reconnector.FinalizeOptimization(removeFreeVertexes)
                        
        refined_vol = self._ComputeVolume(model_part.Elements)
        refined_surface = self._ComputeSurfaceArea(model_part.Conditions)

        print("refined_vol = ",refined_vol)
        self.assertAlmostEqual(refined_vol, original_vol, 12)
        #self.assertAlmostEqual(refined_surface, original_surface, 9)
        self.assertEqual(len(model_part.Elements), 2003)
        self.assertEqual(len(model_part.Nodes), 462)
        print(model_part)


if __name__ == '__main__':
    KratosUnittest.main()
