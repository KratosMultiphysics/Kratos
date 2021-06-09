import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
from matplotlib import pyplot as plt
import numpy as np
import json 
from math import atan2, asin, pi, copysign

# This process reads, modifies and creates vtk files to postprocess the results obtained thorught Simscape Multibody software to complement master stiffnes matrix process
import vtk # Read, modify and create vtk files.
from vtk.numpy_interface import dataset_adapter as dsa # Enables to work with vtk data as python arrays
import numpy as np # This needs to be replaced with Kratos vector and/or matrix
import os # Used to create the folder of results
from scipy.io import loadmat # Load .mat files
import KratosMultiphysics 
from math import atan2, asin, pi, copysign
from multiprocessing import Process
import time 

def QuaternionToEulerAngles(q):
    """
    Convert a quaternion into euler angles (roll, pitch, yaw)
    roll is rotation around x in radians (counterclockwise)
    pitch is rotation around y in radians (counterclockwise)
    yaw is rotation around z in radians (counterclockwise)
    """
    # roll (x-axis rotation)
    sinr_cosp = 2 * (q.W * q.X + q.Y * q.Z)
    cosr_cosp = 1 - 2 * (q.X * q.X + q.Y * q.Y)
    roll = atan2(sinr_cosp, cosr_cosp)

    # pitch (y-axis rotation)
    sinp = 2 * (q.W * q.Y - q.Z * q.X)
    if (abs(sinp) >= 1):
        pitch = copysign(pi / 2, sinp) # use 90 degrees if out of range
    else:
        pitch = asin(sinp)

    # yaw (z-axis rotation)
    siny_cosp = 2 * (q.W * q.Z + q.X * q.Y)
    cosy_cosp = 1 - 2 * (q.Y * q.Y + q.Z * q.Z)
    yaw = atan2(siny_cosp, cosy_cosp)

    angles = []
    angles = np.vstack((roll,pitch,yaw)).T
    return angles[0,0], angles[0,1], angles[0,2]


class RunHROM(StructuralMechanicsAnalysisROM):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)

            ######### GUI attributes
        self.Continue = True
        self.time_step = 0
        self.ReadBCs()
        self.UpdateBC()
        

    def ModifyInitialGeometry(self):
        """Here is the place where the HROM_WEIGHTS are assigned to the selected elements and conditions"""
        super().ModifyInitialGeometry()
        computing_model_part = self._solver.GetComputingModelPart()
        ## Adding the weights to the corresponding elements
        with open('ElementsAndWeights.json') as f:
            HR_data = json.load(f)
            for key in HR_data["Elements"].keys():
                computing_model_part.GetElement(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Elements"][key])
            for key in HR_data["Conditions"].keys():
                computing_model_part.GetCondition(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Conditions"][key])
        
    
    def UpdateBC(self):
        BC = self.getBC()
        hrom_parameters = self.project_parameters["processes"]["constraints_process_list"][0]["Parameters"]
        if not hrom_parameters.Has("boundary_condition"):
            hrom_parameters.AddEmptyArray("boundary_condition")
        hrom_parameters["boundary_condition"].SetVector(BC)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if self.time_step < self.time_steps-1:
            self.time_step += 1
            self.UpdateBC()
        else:
            self.Continue = False
        

    def KeepAdvancingSolutionLoop(self):
        return self.Continue

    
    def getBC(self):
        BC = KratosMultiphysics.Vector(self.DoF)
        for i in range(self.DoF):
            BC[i] = self.BCs[self.time_step,i]
        return BC

    def ReadBCs(self):
        hrom_parameters = self.project_parameters["processes"]["constraints_process_list"][0]["Parameters"]
        self.start_time_step = int(hrom_parameters["time_step_range"].GetVector()[0])
        self.end_time_step = int(hrom_parameters["time_step_range"].GetVector()[1])
        # Read matlab results, now we are proposing them
        # Creating results to validate results
        data = loadmat('sim_res.mat') # Load the simscape results
        Number_of_frames = int(data["Num_of_frames"])-1 # Number of frames -1, we do not consider Frame 1 because it is said to be fixed
        self.DoF = int(Number_of_frames*6) # Degrees of freedom used to post-process
        # self.time_steps = len(data["Time"]) # Time step lenght
        self.time_steps = self.end_time_step-self.start_time_step
        print("Matlab data loaded and readed...")

        self.BCs = KratosMultiphysics.Matrix(self.time_steps,self.DoF) # Matrix of disp/rot relations between Kratos' process and Simulink's multibody analysis
        for i in range(Number_of_frames):
            Interface_name = "Interface_"+str(i+2) # Starts with Interface_2 (All displacements and rotations are measured w.r.t Interface Frame 1)
            pos_x = int(i*6)
            pos_y = pos_x+1
            pos_z = pos_y+1
            pos_quat = int(3*(2*i+1))
            # Saved as dx1,dy1,dz1,rx1,ry1,rz1,dx2,dy2.....rzn
            counter = 0
            for j in range(self.start_time_step,self.end_time_step):
                self.BCs[counter,pos_x] = data[Interface_name]["Disp_x"][0][0][j][0]/1000 # [access to array of arrays]->[access to array]->[time step row]->[access to value]
                self.BCs[counter,pos_y] = data[Interface_name]["Disp_y"][0][0][j][0]/1000
                self.BCs[counter,pos_z] = data[Interface_name]["Disp_z"][0][0][j][0]/1000
                Quat = KratosMultiphysics.Quaternion()
                Quat_data = data[Interface_name]["Quaternion"][0][0][j]
                Quat.W = Quat_data[0]
                Quat.X = Quat_data[1]
                Quat.Y = Quat_data[2]
                Quat.Z = Quat_data[3]
                self.BCs[counter,pos_quat], self.BCs[counter,pos_quat+1], self.BCs[counter,pos_quat+2] = QuaternionToEulerAngles(Quat) # Save euler angles
                counter+=1
        print("Boundary conditions saved...")


##############################################################################################
#                                          RUN HROM                                          #
##############################################################################################
if __name__ == "__main__":
    #KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    with open("ProjectParameters_HROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = RunHROM(model,parameters)
    simulation.Run()
