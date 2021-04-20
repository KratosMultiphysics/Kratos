# This process reads, modifies and creates vtk files to postprocess the results obtained thorught Simscape Multibody software to complement master stiffnes matrix process
import vtk # Read, modify and create vtk files.
from vtk.numpy_interface import dataset_adapter as dsa # Enables to work with vtk data as python arrays
import numpy as np # This needs to be replaced with Kratos vector and/or matrix
import os # Used to create the folder of results
from scipy.io import loadmat # Load .mat files
import KratosMultiphysics 
from multiprocessing import Process
import time 

def CreateVtkFiles(step_in,step_fin,mesh,SnapshotMatrix_stresses,folder_name):
    for i in range(step_in,step_fin):# Loop in time results
        meshNew = dsa.WrapDataObject(mesh) # Create a base mesh (we can take any mesh and reset de values of the displacement)
        von = meshNew.PointData['VON_MISES_STRESS'] # Get the stresses pointer
        for j in range(int(von.shape[0])):
            von[j] = SnapshotMatrix_stresses[j,i]
        # Write a new vtk file every time step with the new results
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(folder_name+"/New"+str(i+1)+".vtk")
        writer.SetInputData(meshNew.VTKObject)
        writer.Write()
        #print("Postprocess for time "+str(i))


# Create new folder for the new vtk output files
def RunPostProcess(index):
    folder_name = 'vtk_postprocess_'+str(index)
    if not os.path.exists(folder_name):
        os.mkdir(folder_name) # Create folder
    else:
        for f in os.listdir(folder_name):
            os.remove(os.path.join(folder_name, f))
    print("Vtk folder created")

    # Read vtk output file
    time_steps = 60
    fileName = 'vtk_output/Structure_0_1.vtk' # Disp/Rot are obtained with respect to Interface Frame 1, first 6 DoF will not be necessary
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    mesh = reader.GetOutput() # Append the mesh in the mesh list

    n = 10
    time_steps_in_open = int(time_steps/n)
    time_res = time_steps%n
    processes = []
    print("Postprocess Started")
    print("Creating Vtk files....")
    start = time.perf_counter()
    with open('new_stress.npy', 'rb') as f:
        SnapshotMatrix_stresses = np.load(f)
    for k in range(n):
        step_in = time_steps_in_open*k
        step_fin = time_steps_in_open*(k+1)
        if k==n-1:
            step_fin += time_res
        p = Process(target=CreateVtkFiles, args=(step_in,step_fin,mesh,SnapshotMatrix_stresses,folder_name))
        processes.append(p)
        p.start()

    for process in processes:
        process.join()
    print("Postprocess finished")
    finish = time.perf_counter()
    print("Finished in ",round(finish-start,2)," second(s)")

