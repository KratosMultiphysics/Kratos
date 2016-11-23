from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *
CheckForPreviousImport()

class MeshAdaptor:
    def __init__(self,model_part):
        self.model_part = model_part
        
    def DoRefinement(self,target_distance, improve_quality=True):
        
        for node in self.model_part.Nodes:
            node.Set(TO_REFINE,False)
            
        for elem in self.model_part.Elements:
            if elem.GetValue(SPLIT_ELEMENT) == True:
                for node in elem.GetNodes():
                    node.Set(TO_REFINE,True)
                        
        for elem in self.model_part.Elements:
            to_split = False
            for node in elem.GetNodes():
                if node.Is(TO_REFINE):
                    to_split = True
                    
            if to_split == True:
                elem.SetValue(SPLIT_ELEMENT,True)
            
        #do refinement
        nodal_neighbour_search = FindNodalNeighboursProcess(self.model_part,20,20)
        nodal_neighbour_search.Execute()

        print("before refinement")
        ###perform the refinement
        Refine = LocalRefineTetrahedraMesh(self.model_part)
        refine_on_reference = False;
        interpolate_internal_variables = False;
        Refine.LocalRefineMesh(refine_on_reference,interpolate_internal_variables)    

        ###recompute the neighbours since they changed due to the creation of new nodes
        nodal_neighbour_search.ClearNeighbours()
        nodal_neighbour_search.Execute()
        print("after refinement")
        ##imporving mesh quality after refining
        #OptimizeQuality(ModelPart& r_model_part,int simIter, int iterations ,
                         #bool processByNode, bool processByFace, bool processByEdge,
                         #bool saveToFile, bool removeFreeVertexes ,
                         #bool evaluateInParallel , bool reinsertNodes , bool debugMode, int minAngle)
                         
        if improve_quality == True:
            reconnector = TetrahedraReconnectUtility(self.model_part)
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
            minAngle = 5


            reconnector.setBlockSize(2048)
            reconnector.OptimizeQuality(self.model_part, simIter, iterations, ProcessByNode, ProcessByFace, ProcessByEdge, saveToFile, removeFreeVertexes, evaluateInParallel, reInsertNodes, debugMode,minAngle)
            meshIsValid = reconnector.EvaluateQuality()
            reconnector.FinalizeOptimization(removeFreeVertexes)
                
