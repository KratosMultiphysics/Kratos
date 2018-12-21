from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
CheckForPreviousImport()

#import matplotlib
import collections

import numpy as np
#from pylab import *

import matplotlib.pyplot as plt

# time control starts
from time import *
print(ctime())
# measure process time
t0p = clock()
# measure wall time
t0w = time()


class InterpretLabTestData:
    #

    def __init__(self, model_part, problem_path):

        # general information of problem
        self.model_part = model_part
        self.problem_path = problem_path
        self.mesh_id = 0
        self.output_files = []
        self.output_files.append(0)
        self.THEarray = []

    # set mesh_id and create output file including header
    def Initialize(self, mesh_id):

        self.mesh_id = mesh_id
        self.elems = self.model_part.GetElements(self.mesh_id)
  
        file_name_out = "LabTestData_" + str(0) + ".csv"
        self.figure_path_out = os.path.join(self.problem_path, file_name_out)

        # write file headers
        if(os.path.exists(self.figure_path_out) == False):
            figure_file = open(self.figure_path_out, "w")
            line_header  = "t ID p q theta p0 b wP epsVol h G e F sigT sigE" +"\n"
            figure_file.write(line_header)
            figure_file.close()

    #
    def SetStepResults(self):

        self._SetStepResultsGauss(self.figure_path_out)

    #
    def _GetStepTime(self):

        return self.model_part.ProcessInfo[TIME]

    #
    def _GetStepDeltaTime(self):

        return self.model_part.ProcessInfo[DELTA_TIME]

    #
    def _SetStepResultsGauss(self, figure_path):

        proc_info = self.model_part.ProcessInfo

        #
        elem = self.elems[1]
            
        # get values for all gauss points of the element
        gaussQuantities = []
        pp = elem.GetValuesOnIntegrationPoints( STRESS_INV_P, proc_info )[0][0]
        qq = elem.GetValuesOnIntegrationPoints( STRESS_INV_J2, proc_info )[0][0]
        p0 = elem.GetValuesOnIntegrationPoints( PRECONSOLIDATION, proc_info )[0][0]
        bb = elem.GetValuesOnIntegrationPoints( BONDING, proc_info )[0][0]
        theta = elem.GetValuesOnIntegrationPoints( STRESS_INV_THETA, proc_info )[0][0]
        # Hencky strain tensor
        FF = self._ConvertTulpeToMatrix(elem.GetValuesOnIntegrationPoints( TOTAL_DEFORMATION_GRADIENT, proc_info )[0])
        hh = self._ComputeHenckyStrainFromF(FF)
        epsVol = hh[0,0] + hh[1,1] + hh[2,2]
        hhSt = self._WriteTulpeToString(hh)
        # stress tensors
        sigT = self._ConvertTulpeToMatrix(elem.GetValuesOnIntegrationPoints( TOTAL_CAUCHY_STRESS, proc_info )[0])
        sigTSt = self._WriteTulpeToString(sigT)
        sigE = self._ConvertTulpeToMatrix(elem.GetValuesOnIntegrationPoints( CAUCHY_STRESS_TENSOR, proc_info )[0])
        sigESt = self._WriteTulpeToString(sigE)
        # water pressure
        ww = self._CalculateWaterPressure(sigT, sigE)
        # ID
        ii = elem.Id
        # strain tensors
        FFSt = self._WriteTulpeToString(FF)
        GG = self._ConvertTulpeToMatrix(elem.GetValuesOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_TENSOR, proc_info )[0])
        GGSt = self._WriteTulpeToString(GG)
        ee = self._ConvertTulpeToMatrix(elem.GetValuesOnIntegrationPoints( ALMANSI_STRAIN_TENSOR, proc_info )[0])
        eeSt = self._WriteTulpeToString(ee)

        # create string
        line_value = str(self._GetStepTime())+ " " + str(ii) + " " + str(pp) + " " + str(qq) + " " + str(theta) + " " + str(p0) + " " + str(bb) + " "
        line_value = line_value + str(ww) + " " + str(epsVol) + " " + hhSt + " " + GGSt + " " + eeSt + " " + FFSt + " " + sigTSt + " " + sigESt + " "
        line_value = line_value + "\n"

        # write line in output file
        figure_file = open(figure_path, "a")
        figure_file.write(line_value)
        figure_file.close()

    #
    def _WriteTulpeToString(self, inputT):

        outString = str(inputT[0,0]) + " " + str(inputT[0,1]) + " " + str(inputT[0,2]) + " "
        outString = outString + str(inputT[1,0]) + " " + str(inputT[1,1]) + " " + str(inputT[1,2]) + " "
        outString = outString + str(inputT[2,0]) + " " + str(inputT[2,1]) + " " + str(inputT[2,2])
        return outString

    #
    def _ConvertTulpeToMatrix(self, inputT):

        FF = Matrix(3,3)
        FF[0,0] = inputT[0]
        FF[0,1] = inputT[1]
        FF[0,2] = inputT[2]
        FF[1,0] = inputT[3]
        FF[1,1] = inputT[4]
        FF[1,2] = inputT[5]
        FF[2,0] = inputT[6]
        FF[2,1] = inputT[7]
        FF[2,2] = inputT[8]
        return FF

    #
    def _CalculateWaterPressure(self, sigT, sigE):

        wwMat = sigT - sigE
        ww = (wwMat[0,0] + wwMat[1,1] + wwMat[2,2])/3
        return ww

    #
    def _ComputeHenckyStrainFromF(self, FF):
        import numpy as np
        #create numpy matrices
        #Fmat = np.matrix(np.identity(3), copy=False)
        Fmat = np.identity(3)
        for i in range (0,3):
            for j in range (0,3):
                Fmat[i,j] = FF[i,j]
        #calculate left Cauchy-Green tensor and its eigenvalues
        bb = np.matmul(Fmat, np.transpose(Fmat))
        ll, vv = np.linalg.eig(bb)
        #modify eigenvalues and write into matrix LL
        LL = np.identity(3)
        for i in range (0,3):
            LL[i,i] = np.log(ll[i])/2.0
        #calculate hencky strain
        hh = np.matmul(vv,LL)
        hh = np.matmul(hh,np.transpose(vv))
        #assemble Kratos matrix and output
        hhMat = Matrix(3,3)
        for i in range (0,3):
            for j in range (0,3):
                hhMat[i,j] = hh[i,j]
        #assemble strain vector in Voigt notation
        hhVec = Vector(6)
        for i in range (0,3):
            hhVec[i] = hh[i,i]
        hhVec[3] = hh[0,1]*2.0
        if(6 == 6):
            hhVec[4] = hh[1,2]*2.0
            hhVec[5] = hh[0,2]*2.0
        return hhMat

    #
    def _ComputeGreenLagrangeFromF(self, FF):
        import numpy as np
        #create numpy matrices
        #Fmat = np.matrix(np.identity(3), copy=False)
        Fmat = np.identity(3)
        for i in range (0,3):
            for j in range (0,3):
                Fmat[i,j] = FF[i,j]
        #calculate right Cauchy-Green tensor
        cc = np.matmul(np.transpose(Fmat), Fmat)
        #calculate Green-Lagrange strain
        gg = cc - np.identity(3)
        gg = gg*0.5
        #assemble Kratos matrix and output
        ggMat = Matrix(3,3)
        for i in range (0,3):
            for j in range (0,3):
                ggMat[i,j] = gg[i,j]
        #assemble strain vector in Voigt notation
        ggVec = Vector(6)
        for i in range (0,3):
            ggVec[i] = gg[i,i]
        ggVec[3] = gg[0,1]*2.0
        if(6 == 6):
            ggVec[4] = gg[1,2]*2.0
            ggVec[5] = gg[0,2]*2.0
        return ggMat

    


            
