import KratosMultiphysics

import sys
import time
import math

import numpy as np
from numpy import linalg as la

class FileWriter:
    def __init__(self, FileName, DataNames, OpenMode="w"):
        if type(DataNames) is not list:
            raise Exception("The result column names have to be passed as list!")
        
        self.file = open(FileName, OpenMode)

        self.num_results = len(DataNames)

        for name in DataNames:
            self.file.write(str(name) + "\t")
        self.file.write("\n")
        self.file.flush()

    def WriteToFile(self, Results):
        if type(Results) is not list:
            raise Exception("The results  have to be passed as list!")
        if self.num_results != len(Results):
            raise Exception("Wrong number of results passed")
        
        for result in Results:
            self.file.write(str(result) + "\t")
        self.file.write("\n")
        self.file.flush()

    def CloseFile(self):
        self._close_file()
        
    def __del__(self): # in case the user forgets to close the file
        self._close_file()

    def _close_file(self):
        try:
            if not self.file.closed:
                self.file.close()
                print("Result file was closed")
        except AttributeError:
            print("File Object did not exist")
        

def TimeRoundValue(DeltaTime):
    return abs(int(math.log10(DeltaTime))) + 2


def Norm(array):
    return la.norm(array)


def SetMeshVelocityToFluid(Nodes, FixDofs=False):
    KratosMultiphysics.VariableUtils().CopyVectorVar(KratosMultiphysics.MESH_VELOCITY, KratosMultiphysics.VELOCITY, Nodes)
    if FixDofs: # This should not be necessary
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_X, True, Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Y, True, Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Z, True, Nodes)

    # OLD:
    # for node in Nodes:
    #     mesh_velocity = node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, 0)  # get mesh velocity at current step
    #     node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, mesh_velocity)  # assign it to fluid velocity at current step
    #     node.Fix(KratosMultiphysics.VELOCITY_X)
    #     node.Fix(KratosMultiphysics.VELOCITY_Y)
    #     node.Fix(KratosMultiphysics.VELOCITY_Z)
    

def GetDisplacements(NodesOfStructure, Dimension=3):
    #displacements = [0.0]*3*len(NodesOfStructure)
    displacements = np.zeros(3 * len(NodesOfStructure)) 
    index = 0
    for node in NodesOfStructure:
        displacements[3*index] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0)
        displacements[3*index+1] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0)
        if Dimension == 3:
            displacements[3*index+2] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
        index += 1


    return displacements

def SetDisplacements(displacements, NodesOfStructure, Dimension=3):
    index = 0
    for node in NodesOfStructure:
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0,displacements[3*index])
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,displacements[3*index+1])
        if Dimension == 3:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0,displacements[3*index+2])
            
        index += 1

def CalculateResidual(Solution, Old_Solution):
    #residual = [0.0]*len(Solution)
    residual = np.zeros(len(Solution))
    for index in range(len(Solution)):
        residual[index] = Solution[index] - Old_Solution[index]
    return residual

def CalculateRelaxation(RelaxationCoefficient, Old_Solution, Residual):
    for index in range(len(Old_Solution)):
        Old_Solution[index] = Old_Solution[index] + RelaxationCoefficient * Residual[index]
    return Old_Solution # this is the relaxed new solution

def ComputeAitkenRelaxation(OldCoefficient, residual, old_residual, iteration):
    # reference for implementation see header
    MaxInitialCoefficient = 0.125 # maximum relaxation coefficient for first iteration
    #print("k = ", iteration)
    if iteration < 1:
        NewCoefficient = min(OldCoefficient,MaxInitialCoefficient)
    else:
        numerator = 0
        denominator = 0
        for i in range(len(residual)):
            numerator += old_residual[i] * (residual[i] - old_residual[i])
            denominator += pow(residual[i]-old_residual[i],2)
        NewCoefficient = - OldCoefficient * (numerator/denominator)
        #if NewCoefficient > 20:
        #    NewCoefficient = 20
        #elif NewCoefficient <= 0:
        #    NewCoefficient = 0.001
    return NewCoefficient






def ApplyVelocityRampUp(model_part, time):
    v_max = 0.06067
    v_t = v_max * 0.50 * (1 - math.cos(3.41 * time *0.10))
    for node in model_part.Nodes:
        h = 0.50
        velocity = -v_t * math.pow(node.Y,2) / math.pow(h,2) + 2 * v_t * node.Y / h
        node.SetSolutionStepValue(VELOCITY_X,0,velocity)
def ApplyVelocityMaximum(model_part):
    v_max = 0.06067
    for node in model_part.Nodes:
        velocity = v_max
        node.SetSolutionStepValue(VELOCITY_X,0,velocity)
