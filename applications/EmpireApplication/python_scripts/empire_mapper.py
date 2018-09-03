from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
import ctypes as ctp
import os

CheckForPreviousImport()

# Import the mapper library
# A header file is put in the header file folder to give the users the syntax for the mapping routines.
# TODO Aditya does the first case work all the time? Also should we change the two cases?
try: # OpenMPI
    libMapper = ctp.CDLL(os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE'], ctp.RTLD_GLOBAL)
    print("::EMPIRE:: Using standard OpenMPI")
except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
    libMapper = ctp.cdll.LoadLibrary(os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE'])
    print("::EMPIRE:: Using Intel MPI or OpenMPI compiled with \"–disable-dlopen\" option")

## Wrapper class for the mapper
# Consturctor will have two model parts as arguments, type of mapper and then options to the mapper. 
# These options should be explained in detail here.
class FiniteElementMapper:
    def __init__(self, name, dimension,  modelPart1, modelPart2, dual=False, oppositeSurfaceNormal=False, enforceConsistency=True):
        self.modelPart1 = modelPart1.Name
        self.modelPart2 = modelPart2.Name
        self.name = name
        self.dimension = dimension

        self.meshModelPart1 = ModelPart(self.modelPart1)
        self.meshModelPart2 = ModelPart(self.modelPart2)

        self.wrapProc1  = WrapperProcess(modelPart1, self.meshModelPart1, self.dimension)
        self.wrapProc2  = WrapperProcess(modelPart2, self.meshModelPart2, self.dimension)

        self.wrapProc1.ExtractInterface()
        self.wrapProc2.ExtractInterface()

       
        numNodes1 = [];         numNodes2 = []
        numElems1 = [];         numElems2 = []
        nodes1     = [];        nodes2    = []
        nodeIDs1   = [];        nodeIDs2  = []
        numNodesPerElem1 = [];  numNodesPerElem2 = []
        elems1         = [];    elems2    = []


        self.wrapProc1.ExtractMeshInfo(numNodes1,numElems1,nodes1,nodeIDs1,numNodesPerElem1,elems1)
        self.wrapProc2.ExtractMeshInfo(numNodes2,numElems2,nodes2,nodeIDs2,numNodesPerElem2,elems2)
        
        print("Extracted Mesh 1 ::\n")
        print("Number Of Nodes    :: ",numNodes1)
        print("Number Of Elements :: ",numElems1)
        #print("Nodes              :: ",nodes1)
        #print("NodeIDs1           :: ",nodeIDs1)
        #print("Number of Nodes Per Elem :: ",numNodesPerElem1)
        #print("Elements           :: ", elems1)
        print("Dimension          :: ", self.dimension)
        
        print("Extracted Mesh 2 ::\n")
        print("Number Of Nodes    :: ",numNodes2)
        print("Number Of Elements :: ",numElems2)
        #print("Nodes              :: ",nodes2)
        #print("NodeIDs1           :: ",nodeIDs2)
        #print("Number of Nodes Per Elem :: ",numNodesPerElem2)
        #print("Elements           :: ", elems2)

        c_numNodes1 = (ctp.c_int * len(numNodes1))(*numNodes1)
        c_numElems1 = (ctp.c_int * len(numElems1))(*numElems1)
        c_nodes1    = (ctp.c_double * len(nodes1))(*nodes1)
        c_nodeIDs1  = (ctp.c_int * len(nodeIDs1)) (*nodeIDs1)
        c_numNodesPerElem1 = (ctp.c_int * len(numNodesPerElem1)) (*numNodesPerElem1)
        c_elems1     = (ctp.c_int * len(elems1)) (*elems1)
        
        c_numNodes2 = (ctp.c_int * len(numNodes2))(*numNodes2)
        c_numElems2 = (ctp.c_int * len(numElems2))(*numElems2)
        c_nodes2    = (ctp.c_double * len(nodes2))(*nodes2)
        c_nodeIDs2  = (ctp.c_int * len(nodeIDs2)) (*nodeIDs2)
        c_numNodesPerElem2 = (ctp.c_int * len(numNodesPerElem2)) (*numNodesPerElem2)
        c_elems2     = (ctp.c_int * len(elems2)) (*elems2)

        c_oppositeSurfaceNormal = (ctp.c_int * 1) (oppositeSurfaceNormal)
        c_enforceConsistency    = (ctp.c_int * 1) (enforceConsistency)
        c_dual    = (ctp.c_int * 1) (dual)

        libMapper.init_FE_MortarMapper(self.name.encode('utf-8'), c_numNodes1[0], c_numElems1[0], c_numNodesPerElem1, c_nodes1, c_nodeIDs1, c_elems1, c_numNodes2[0], c_numElems2[0], c_numNodesPerElem2, c_nodes2, c_nodeIDs2, c_elems2, c_oppositeSurfaceNormal[0], c_dual[0], c_enforceConsistency[0])


# Then we have a do mapping function which takes the two model parts and their fields names as defined in KRATOS as arguments. 
# This function will check in which order the mapping has to be done. 
# No model parts are stored in this class. only names are stored. Then comparison to the names is done. 
    def doMapping(self, modelPart1, dataField1, modelPart2, dataField2):
        df1 = []
        df2 = []

        for node in modelPart1.Nodes:
            val = node.GetSolutionStepValue(dataField1)
            df1.append(val[0])
            df1.append(val[1])
            df1.append(val[2])
           
        for node in modelPart2.Nodes:
            val = node.GetSolutionStepValue(dataField2)
            df2.append(val[0])
            df2.append(val[1])
            df2.append(val[2])
              
        lenDf1 = []
        lenDf2 = []
        lenDf1.append(len(df1))
        lenDf2.append(len(df2)) 

        c_df1 = (ctp.c_double *len(df1)) (*df1)
        c_df2 = (ctp.c_double *len(df2)) (*df2)
        c_lenDf1 = (ctp.c_int * 1 )  (*lenDf1)
        c_lenDf2 = (ctp.c_int * 1 )  (*lenDf2)
        #print(c_df1[0])
        #print(c_df1[1])
        #print(c_df1[2])

        if(modelPart1.Name == self.modelPart1 and modelPart2.Name == self.modelPart2):
            libMapper.doConsistentMapping(self.name.encode('utf-8'), self.dimension, c_lenDf1[0], c_df1, c_lenDf2[0], c_df2)
            i=0
            val = Vector(3)
            for node in modelPart2.Nodes:
                val[0] = c_df2[3*i + 0]
                val[1] = c_df2[3*i + 1]
                val[2] = c_df2[3*i + 2]
                node.SetSolutionStepValue(dataField2,0,val)
                i = i + 1

                     
        if(modelPart1.Name == self.modelPart2 and modelPart2.Name == self.modelPart1):
            libMapper.doConservativeMapping(self.name.encode('utf-8'), self.dimension, c_lenDf1[0], c_df1, c_lenDf2[0], c_df2)
            i=0
            val = Vector(3)
            #print(len(modelPart2.Nodes))
            #print(len(c_df1))
            #print(len(c_df2))
            for node in modelPart2.Nodes:
                val[0] = c_df2[3*i + 0]
                val[1] = c_df2[3*i + 1]
                val[2] = c_df2[3*i + 2]
                node.SetSolutionStepValue(dataField2,0,val)
                i = i + 1
        #print(c_df2[0])
        #print(c_df2[1])
        #print(c_df2[2])         
 
