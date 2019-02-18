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


class InterpretInterestingData:
    #

    def __init__(self, model_part, problem_path):

        # general information of problem
        self.model_part = model_part
        self.problem_path = problem_path
        self.mesh_id = 0
        self.output_files = []
        self.output_files.append(2.5)
        self.output_files.append(2.75)
        self.output_files.append(2.8)
        self.output_files.append(3.05)
        self.output_files.append(3.1)
        self.output_files.append(7.8)
        self.output_files.append(8.05)
        self.output_files.append(8.1)
        self.output_files.append(12.8)
        self.output_files.append(13.05)
        self.output_files.append(13.1)
        self.output_files.append(0.2)
        self.output_files.append(0.35)
        self.output_files.append(0.005)
        self.THEarray = []

    # set mesh_id and create output file including header
    def Initialize(self, mesh_id):

        self.mesh_id = mesh_id
  
        for jj in range(0, len(self.output_files)):
            outputId = self.output_files[jj]
            file_name_nodes = "InterestingNodalData_t" + str(jj) + ".csv"
            file_name_gauss = "InterestingGaussData_t" + str(jj) + ".csv"
            figure_path_nodes = os.path.join(self.problem_path, file_name_nodes)
            figure_path_gauss = os.path.join(self.problem_path, file_name_gauss)

            # write file headers
            if(os.path.exists(figure_path_nodes) == False):
                figure_file = open(figure_path_nodes, "w")
                line_header  = "ID X Y Z disp_X disp_Y disp_Z J WP dWP" +"\n"
                figure_file.write(line_header)
                figure_file.close()

            if(os.path.exists(figure_path_gauss) == False):
                figure_file = open(figure_path_gauss, "w")
                line_header  = "ID X Y Z n1ID n1X n1Y n1Z n2ID n2X n2Y n2Z n3ID n3X n3Y n3Z P J2 Theta Exx Exy Eyx Eyy exx exy exz eyx eyy eyz ezx ezy ezz Fxx Fxy Fxz Fyx Fyy Fyz Fzx Fzy Fzz Qx Qy Qz" +"\n"
                figure_file.write(line_header)
                figure_file.close()

            # write relevant data into THEarray
            auxArray = [self.output_files[jj], figure_path_nodes, figure_path_gauss]
            self.THEarray.append(auxArray)


    #
    def GetStepTime(self):

        return self.model_part.ProcessInfo[TIME]

    #
    def GetStepDeltaTime(self):

        return self.model_part.ProcessInfo[DELTA_TIME]

    #
    def CheckTime(self):

        eps = self.GetStepDeltaTime()/100
        lowerBoundTime = self.GetStepTime() - eps
        upperBoundTime = self.GetStepTime() + eps

        checkArray = [0, 0, 0, False]
        for kk in range(0, len(self.THEarray)):
            if (self.THEarray[kk][0] > lowerBoundTime) and (self.THEarray[kk][0] < upperBoundTime):
                checkTime = True
                checkArray = self.THEarray[kk]
                checkArray.append(checkTime)

        return checkArray

    #
    def SetStepResultsNodes(self, figure_path):

        # start time measurement
        time_start = clock()

        # loop over nodes
        for node in self.model_part.GetNodes(self.mesh_id):
            # get nodal values
            nodalQuantities = []
            nodalQuantities.append(node.Id)
            nodalQuantities.append(node.X)
            nodalQuantities.append(node.Y)
            nodalQuantities.append(node.Z)
            nodalQuantities.append(node.GetSolutionStepValue( DISPLACEMENT_X ))
            nodalQuantities.append(node.GetSolutionStepValue( DISPLACEMENT_Y ))
            nodalQuantities.append(node.GetSolutionStepValue( DISPLACEMENT_Z ))
            nodalQuantities.append(node.GetSolutionStepValue( JACOBIAN ))
            nodalQuantities.append(node.GetSolutionStepValue( WATER_PRESSURE ))
            nodalQuantities.append(node.GetSolutionStepValue( EXCESS_WATER_PRESSURE ))

            # create string
            line_value = ""
            for ii in range(0, len(nodalQuantities)):
                line_value = line_value + str(nodalQuantities[ii]) + " "
            line_value = line_value + "\n"

            # write line in output file
            figure_file = open(figure_path, "a")
            figure_file.write(line_value)
            figure_file.close()

        # end time measurement
        time_end = clock()
        used_time = time_end - time_start
        line_value = "time step = " + str(self.GetStepTime()) + "s  --  process time = " + str(round(used_time,2)) + "s" + "\n"
        figure_file = open(figure_path, "a")
        figure_file.write(line_value)
        figure_file.close()

    #
    def SetStepResultsGauss(self, figure_path):

        # start time measurement
        time_start = clock()

        proc_info = self.model_part.ProcessInfo

        # loop over elements
        for elem in self.model_part.GetElements(self.mesh_id):
            elementNodes = []
            # loop over element nodes
            for node in elem.GetNodes():
                elementStuff = []
                elementStuff.append(node.Id)
                elementStuff.append(node.X)
                elementStuff.append(node.Y)
                elementStuff.append(node.Z)
                elementNodes.append(elementStuff)
                
            # get values for all gauss points of the element
            gaussQuantities = []
            gaussQuantities.append(elem.Id)
            gaussQuantities.append(elem.GetIntegrationPoints())
            gaussQuantities.append(elementNodes)
            gaussQuantities.append(elem.GetValuesOnIntegrationPoints( STRESS_INV_P, proc_info ))
            gaussQuantities.append(elem.GetValuesOnIntegrationPoints( STRESS_INV_J2, proc_info ))
            gaussQuantities.append(elem.GetValuesOnIntegrationPoints( STRESS_INV_THETA, proc_info ))
            gaussQuantities.append(elem.GetValuesOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_TENSOR, proc_info ))
            gaussQuantities.append(elem.GetValuesOnIntegrationPoints( ALMANSI_STRAIN_TENSOR, proc_info ))
            gaussQuantities.append(elem.GetValuesOnIntegrationPoints( TOTAL_DEFORMATION_GRADIENT, proc_info ))
            #gaussQuantities.append(elem.GetValuesOnIntegrationPoints( DEFORMATION_GRADIENT, proc_info ))
            gaussQuantities.append(elem.GetValuesOnIntegrationPoints( DARCY_FLOW, proc_info ))

            # create string
            line_value = ""
            for ii in range(0, len(gaussQuantities)):
                line_value = line_value + str(gaussQuantities[ii]) + " "
            line_value = line_value + "\n"

            # write line in output file
            figure_file = open(figure_path, "a")
            figure_file.write(line_value)
            figure_file.close()

        # end time measurement
        time_end = clock()
        used_time = time_end - time_start
        line_value = "time step = " + str(self.GetStepTime()) + "s  --  process time = " + str(round(used_time,2)) + "s" + "\n"
        figure_file = open(figure_path, "a")
        figure_file.write(line_value)
        figure_file.close()

    #
    def SetStepResults(self):

        checkArray = self.CheckTime()
        if ( checkArray[3] == True):
            self.SetStepResultsNodes(os.path.join(self.problem_path, checkArray[1]))
            self.SetStepResultsGauss(os.path.join(self.problem_path, checkArray[2]))


            
