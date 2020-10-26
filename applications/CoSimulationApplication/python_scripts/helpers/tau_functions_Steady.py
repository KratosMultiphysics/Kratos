# -*- coding: utf-8 -*-
# make script backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import re, glob, subprocess, time, os, warnings, json
import numpy as np
try:
    import tau_python
    from tau_python import tau_msg
    from tau_python import tau_solver_unsteady_get_physical_time
    from tau_python import tau_mpi_rank
    rank = tau_mpi_rank()
    from tau_python import *
    import PyPara, PySurfDeflect
    tau_available = True
except:
    tau_available = False
    warnings.warn('tau modules not available')

from scipy.io import netcdf
import ChangeFormat as CF

# Assign default settings
working_path = os.getcwd() + '/'
#path = "/work/piquee/MembraneWing/kratos_fsi_big"
with open(working_path + 'input/tau_settings.json') as json_file:
    tau_settings = json.load(json_file)

start_step = tau_settings["start_step"]
echo_level = tau_settings["echo_level"]
rotate = tau_settings["rotate"]

# Read settings
try:
    with open(working_path + '/input/tau_settings.json') as json_file:
        tau_settings = json.load(json_file)

    start_step = tau_settings["start_step"]
    tau_path = tau_settings["tau_path"]
    echo_level = tau_settings["echo_level"]
    rotate = tau_settings["rotate"]
except:
    warnings.warn('/input/tau_settings.json not found')



# Remove output files and deform mesh files from previous simulations
def RemoveFilesFromPreviousSimulations():
    # Get a list of all the output and mesh file paths
    file_list = glob.glob(working_path + 'Outputs/*')
    file_list += glob.glob(working_path + 'Mesh/airfoil_Structured_scaliert.grid.def*')

    # Iterate over the list of filepaths and remove each file.
    for file_path in file_list:
        try:
            os.remove(file_path)
        except:
            raise Exception('Error while deleting file : "{}" '.format(file_path))


# Convert tau output to dat file using tau2plt
def ConvertOutputToDat(working_path, tau_path, step, para_path_mod, ouput_file_pattern, step_mesh):
    PrintBlockHeader("Start Garthering Solution Data at time %s" % (str(time)))
    # Execute gather
    command = tau_path + 'gather ' +  para_path_mod
    print(command)
    subprocess.call(command, shell=True)

    PrintBlockHeader("Start Writting Solution Data at time %s" % (str(time)))
    # Write Tautoplt.cntl file
    tautoplt_filename = WriteTautoplt(working_path, step, para_path_mod, ouput_file_pattern, step_mesh)

    # Execute tau2plt to convert output file into dat
    command = tau_path + 'tau2plt ' + tautoplt_filename
    subprocess.call(command, shell=True)
    PrintBlockHeader("Stop Writting Solution Data at time %s" % (str(time)))


def PrintBlockHeader(header):
    tau_python.tau_msg("\n" + 50 * "*" + "\n" + "* %s\n" %header + 50*"*" + "\n")


# Write membrane's displacments in a file
def ExecuteBeforeMeshDeformation(total_displacements, step, para_path_mod, start_step):
    # Compute relative displacements
    relative_displacements = ComputeRelativeDisplacements(total_displacements, step, start_step)
    return relative_displacements
    # Read tau's parameter file
#    Para = PyPara.Parafile(para_path_mod)

    # Read the interface fluid grid
#    print("before PySurf")
#    ids, coordinates = PySurfDeflect.read_tau_grid(Para)
#    print("after PySurf")
    # Write membrane's displacments in a file
    #if tau_mpi_rank() == 0:
#    WriteInterfaceDeformationFile(ids, coordinates, relative_displacements)


# Computes fluid forces at the nodes
def ComputeFluidForces(working_path, step, word, ouput_file_pattern):
    # Read mesh and pressure from interface file
    X, Y, Z, nodal_pressures, elem_connectivities = ReadTauOutput(working_path, step, 41.09, word, ouput_file_pattern)

    # calculating the force vector
    fluid_forces = CalculateNodalFluidForces(X, Y, Z, nodal_pressures, elem_connectivities)

    return fluid_forces


# GetFluidMesh is called only once at the beginning, after the first fluid solve
def GetFluidMesh(working_path, step, word, ouput_file_pattern):
    # Read mesh from interface file
    X, Y, Z, P, elem_connectivities = ReadTauOutput(working_path, step, 41.09, word, ouput_file_pattern)
    print("after ReadTauOutput")

    # Transform nodal coordinates to numpy array
    nodal_coords = ReadNodalCoordinates(X, Y, Z)
    print("after ReadNodalCoordinates")
    # Save element types in a numpy array
    element_types = ReadElementTypes(int(len(elem_connectivities)/4))

    # In vtk format element connectivities start from 0, not from 1
    elem_connectivities -= 1

    return nodal_coords, elem_connectivities, element_types


# Write Tautoplt.cntl file
def WriteTautoplt(working_path, step, para_path_mod, ouput_file_pattern, step_mesh):
    # Define Tautoplt.cntl file name and check if it already exists
    tautoplt_filename = working_path + 'Tautoplt.cntl'
    RemoveFileIfExists(tautoplt_filename)
    initial_tautoplt_filename = working_path + 'Tautoplt_initial.cntl'
    CheckIfPathExists(initial_tautoplt_filename)

    # Read and write simultaneously
    tautoplt_file_writing = open(tautoplt_filename, 'w')
    with open(initial_tautoplt_filename, 'r+') as tautoplt_file_reading:
        # Loop over lines
        for line in tautoplt_file_reading:
            # Check if line is from IO section and modify it
            line = ModifyFilesIOLines(line, working_path, step, para_path_mod, ouput_file_pattern, step_mesh)
            tautoplt_file_writing.write(line)

        # Close files
        tautoplt_file_writing.close()
        tautoplt_file_reading.close()

    return tautoplt_filename


# Compute relative displacements
def ComputeRelativeDisplacements(total_displacements, step, start_step):
    global previous_total_displacements

    if(step == start_step):
        previous_total_displacements = np.zeros(len(total_displacements))

    # Declaring and initializing relative_displacements vector
    number_of_nodes = int(len(total_displacements)/3)
    relative_displacements = np.zeros([number_of_nodes, 3])

    # Loop over nodes
    for node in range(number_of_nodes):
        # Loop over xyz components
        for j in range(3):
            # Compute relative displacement
            relative_displacements[node, j] = total_displacements[3*node+j]
            relative_displacements[node, j] -= previous_total_displacements[3*node+j]

    previous_total_displacements = total_displacements

    return relative_displacements


# Change format displacements
def ChangeFormatDisplacements(total_displacements):

    # Declaring and initializing new displacement vector
    number_of_nodes = int(len(total_displacements)/3)
    new_displacements = np.zeros([number_of_nodes, 3])

    # Loop over nodes
    for node in range(number_of_nodes):
        # Loop over xyz components
        for j in range(3):
            # Compute relative displacement
            new_displacements[node, j] = total_displacements[3*node+j]

    return new_displacements

# Write membrane's displacments in a file
def WriteInterfaceDeformationFile(ids, coordinates, relative_displacements, word):
    # Open interface_deformfile
    ncf = netcdf.netcdf_file(word + '.interface_deformfile.nc', 'w')

    # define dimensions
    nops = 'no_of_points'
    number_of_points = len(ids[:])
    print('len(ids[:] = ', str(len(ids[:])))
    print('len(relative_displacements) = ', str(len(relative_displacements)))
    ncf.createDimension(nops, number_of_points)

    # define variables
    global_node_ids = ncf.createVariable('global_id', 'i', (nops,))
    nodal_x_coordinates = ncf.createVariable('x', 'd', (nops,))
    nodal_y_coordinates = ncf.createVariable('y', 'd', (nops,))
    nodal_z_coordinates = ncf.createVariable('z', 'd', (nops,))
    nodal_x_displacements = ncf.createVariable('dx', 'd', (nops,))
    nodal_y_displacements = ncf.createVariable('dy', 'd', (nops,))
    nodel_z_displacements = ncf.createVariable('dz', 'd', (nops,))

    # write data
    global_node_ids[:] = ids
    nodal_x_coordinates[:] = coordinates[0,:]
    nodal_y_coordinates[:] = coordinates[1,:]
    nodal_z_coordinates[:] = coordinates[2,:]
    nodal_x_displacements[:] = relative_displacements[:,0]
    nodal_y_displacements[:] = relative_displacements[:,1]
    nodel_z_displacements[:] = relative_displacements[:,2]
    ncf.close()

def ChangeFormat(working_path, step, word1, word2, ouput_file_pattern):
    # Find the interface file name
    interface_filename = FindInterfaceFilename(working_path, step, ouput_file_pattern)
    # Change Format of 'interface_filename'
    print("Looking for number 1")
    number1 = CF.Countline(interface_filename, word1)
    print("number 1 = ", str(number1))
    print("Looking for number 1")
    number2 = CF.Countline(interface_filename, word2)
    print("number 2 = ", str(number2))
    print("Looking for endFile")
    endline = CF.Endline(interface_filename)
    print("endline = ", str(endline))

    with open(interface_filename,'r') as fread:
        with open(interface_filename[0:len(interface_filename)-4] + '.' + word1 + '.dat','w') as fwrite_word1:
            with open(interface_filename[0:len(interface_filename)-4] + '.' + word2 + '.dat','w') as fwrite_word2:
                header1 = fread.readline()
                fwrite_word1.write(header1)
                fwrite_word2.write(header1)
                header2 = fread.readline()
                fwrite_word1.write(header2)
                fwrite_word2.write(header2)

    if number1 < number2:
        CF.WriteNewFile(interface_filename, interface_filename[0:len(interface_filename)-4] + '.' + word1 + '.dat', number1, number2)
        CF.WriteNewFile(interface_filename, interface_filename[0:len(interface_filename)-4] + '.' + word2 + '.dat', number2, endline)
    else:
        CF.WriteNewFile(interface_filename, interface_filename[0:len(interface_filename)-4] + '.' + word2 + '.dat', number2, number1)
        CF.WriteNewFile(interface_filename, interface_filename[0:len(interface_filename)-4] + '.' + word1 + '.dat', number1, endline)


# Read mesh and data from tau output file
def ReadTauOutput(working_path, step, velocity, word, ouput_file_pattern):
    # Find the interface file name
    interface_filename_original = FindInterfaceFilename(working_path, step, ouput_file_pattern)
    print("interface_filename_original = ",interface_filename_original)
    interface_filename = interface_filename_original[0:len(interface_filename_original)-4] + '.' + word + '.dat'

    # Read interface file
    position_info, mesh_info, nodal_data, elem_connectivities = ReadInterfaceFile(
        interface_filename)
    print("after ReadInterfaceFile")
    # Read mesh info
    NodesNr = mesh_info[0]
    print("NodesNr = ", NodesNr)
    X, Y, Z = SaveCoordinatesList(nodal_data, position_info, NodesNr)
    print("SaveCoordinatesList")
    P = SavePressure(nodal_data, position_info, NodesNr, velocity)
    print("SavePressure")

    return X, Y, Z, P, elem_connectivities


# Calculate the fluid forces at the nodes
def CalculateNodalFluidForces(X, Y, Z, nodal_pressures, elem_connectivities):
    nodal_forces = np.zeros(3*len(X))
    # Loop over cells
    for cell in range(int(len(elem_connectivities)/4)):
        # Get the node ids of the cell
        node_ids = GetCellNodeIds(elem_connectivities, cell)

        # Calculate cell force
        cell_force = CalculateCellForce(node_ids, nodal_pressures, X, Y, Z)

        # Extrapolating cell force to the nodes
        for node in range(4):
            # Loop over xyz components
            for component in range(3):
                nodal_forces[3*node_ids[node]+component] += 0.25 * cell_force[component]

    return nodal_forces


# Transform nodal coordinates to numpy array
def ReadNodalCoordinates(X, Y, Z):
    # array to store the coordinates of the nodes: x1,y1,z1,x2,y2,z2,...
    nodal_coords = np.zeros(3*len(X))

    # Loop over nodes
    for node in range(len(X)):
        nodal_coords[3*node+0] = X[node]
        nodal_coords[3*node+1] = Y[node]
        nodal_coords[3*node+2] = Z[node]

    return nodal_coords


# Save element types in a numpy array
def ReadElementTypes(ElemsNr):
    return np.full(ElemsNr, 9, dtype=int)


# Check if file exist and remove it, otherwise print a warning
def RemoveFileIfExists(path):
    if os.path.exists(path):
        os.remove(path)
    else:
        msg = 'The file ' + path + ' does not exist.'
        warnings.warn(msg)


# Check if the path exists
def CheckIfPathExists(path):
    if not os.path.exists(path):
        raise Exception('Path: "{}" not found'.format(path))


# Checks if the line is from IO section in tau2plot.cntl and modifies it
def ModifyFilesIOLines(line, working_path, step, para_path_mod, ouput_file_pattern, step_mesh):
    if 'Primary grid filename:' in line:
        primary_grid_filename = FindPrimaryGridFilename(working_path, step_mesh)
        line = 'Primary grid filename:' + primary_grid_filename + ' \n'
    elif 'Boundary mapping filename:' in line:
        parameter_filename = working_path + para_path_mod
        line = 'Boundary mapping filename:' + parameter_filename + ' \n'
    elif 'Restart-data prefix:' in line:
        output_filename = FindOutputFilename(working_path, step, ouput_file_pattern)
        CheckIfPathExists(output_filename)
        line = 'Restart-data prefix:' + output_filename + ' \n'

    return line


# Find the interface file name
def FindInterfaceFilename(working_path, step, ouput_file_pattern):
    output_filename = FindOutputFilename(working_path, step, ouput_file_pattern)
    interface_filename = output_filename.replace('pval', 'surface.pval')
    if '.dat' not in interface_filename:
        interface_filename = interface_filename.replace('+', '') + '.dat'
    CheckIfPathExists(interface_filename)
    return interface_filename


# Read interface file
def ReadInterfaceFile(interface_filename):
    with open(interface_filename, 'r') as interface_file:
        line = interface_file.readline()
        position_info, mesh_info, line = ReadHeader(interface_file, line)
        nodal_data, line = ReadNodalData(interface_file, line)
        elem_connectivities = ReadElementConnectivities(interface_file, line, ElemsNr=mesh_info[1])
    return position_info, mesh_info, nodal_data, elem_connectivities


# Saves the coordinates into lists
def SaveCoordinatesList(nodal_data, position_info, NodesNr):
    # Read coordinates positions and save them in a list
    x_position = position_info.index('"x"') - 2
    X = nodal_data[(x_position)*NodesNr:(x_position+1)*NodesNr]
    y_position = position_info.index('"y"') - 2
    Y = nodal_data[(y_position)*NodesNr:(y_position+1)*NodesNr]
    z_position = position_info.index('"z"') - 2
    Z = nodal_data[(z_position)*NodesNr:(z_position+1)*NodesNr]
    return X, Y, Z

# Saves the nodal pressure into a numpy array
def SavePressure(nodal_data, position_info, NodesNr, velocity):
    # Read pressure coefficient
    cp_position = position_info.index('"cp"') - 2
    CP = nodal_data[(cp_position)*NodesNr:(cp_position+1)*NodesNr]
    # Compute pressure from cp
    P = [x*velocity*velocity*0.5*1.225 for x in CP]
    P = np.array(P)
    return P


# Get the node ids of the cell
def GetCellNodeIds(elem_connectivities, cell):
    node_ids = np.zeros(4, dtype=int)
    # Loop over cell nodes
    for node in range(4):
        node_ids[node] = elem_connectivities[cell*4+node]-1

    return node_ids


# Calculate cell force
def CalculateCellForce(node_ids, nodal_pressures, X, Y, Z):
    # Calculate cell pressure
    pressure = CalculateCellPressure(nodal_pressures, node_ids)

    # Calculate cell area and normal
    area = CalculateCellArea(X, Y, Z, node_ids)
    normal = CalculateCellNormal(X, Y, Z, node_ids)

    # Calculate cell force
    cell_force = pressure * area * normal

    return cell_force
def CalculateCellArea(X, Y, Z, node_ids):
    # Calculate cell sides
    cell_side_01 = CalculateDistanceVector(X, Y, Z, node_ids[0], node_ids[1])
    cell_side_03 = CalculateDistanceVector(X, Y, Z, node_ids[0], node_ids[3])

    # Calculate cell area
    cell_area = np.cross(cell_side_01, cell_side_03)

    # hold only for 2d case, MUST be modified for 3d interfaces
    cell_area = np.linalg.norm(cell_area)

    return cell_area


# Calculate cell normal
def CalculateCellNormal(X, Y, Z, node_ids):
    # Calculate cell diagonals
    cell_diagonal_02 = CalculateDistanceVector(X,Y,Z,node_ids[0],node_ids[2])
    cell_diagonal_13 = CalculateDistanceVector(X,Y,Z,node_ids[1],node_ids[3])

    # Calculate cell normal
    cell_normal = np.cross(cell_diagonal_02, cell_diagonal_13)

    # Normalize to make unit normal
    magnitude = np.linalg.norm(cell_normal)
    unit_cell_normal = cell_normal/magnitude

    return unit_cell_normal

# Finds steps' corresponding primary grid filename
def FindPrimaryGridFilename(working_path, step_mesh):
    mesh_path = working_path + "Mesh/"
    CheckIfPathExists(mesh_path)
    print('step_mesh = ', step_mesh)
    if step_mesh == 0:
        pattern = 'airfoil_Structured_scaliert.grid'
        return FindInitialMeshFilename(mesh_path, pattern)
    else:
        pattern = 'airfoil_Structured_scaliert.grid.def.'
        mesh_file = FindMeshFilename(mesh_path, pattern, step_mesh)
        if 'domain' in mesh_file:
            position = mesh_file.find('_domain_')
            mesh_filename = mesh_file[0:position]
        else:
            mesh_filename = mesh_file
        print(mesh_filename)
        return mesh_filename


# Finds output filename
def FindOutputFilename(working_path, step, ouput_file_pattern):
    outputs_path = working_path + "Outputs/"
    CheckIfPathExists(outputs_path)
    # if rotate and unsteady "airfoilSol.pval.deform_i="
    # if non rotate and unsteady "airfoilSol.pval.unsteady_i="
    # if rotate and steady
    print("step_FindOutputFilename = ", step)
    CheckIfPathExists(FindFilename(outputs_path, ouput_file_pattern, step))
    return FindFilename(outputs_path, ouput_file_pattern, step)


# Read header from interface file
def ReadHeader(interface_file,line):
    # reading the five lines of the header
    for _ in range(4):
        # reading the position of the coordinates and cp in the list
        if "\"x\" \"y\" \"z\"" in line:
            position_info = line.split()
        # reading the number of nodes and elements
        elif "N=" in line:
            mesh_info = [int(integers)
                         for integers in re.findall(r'\b\d+\b', line)]
        line = interface_file.readline()

    return position_info, mesh_info, line


# Read nodal data from interface file
def ReadNodalData(interface_file,line):
    # reading nodal data
    nodal_data = []
    while 'E+' in line or 'E-' in line:
        for elem in line.split():
            nodal_data.append(float(elem))
        line = interface_file.readline()

    return nodal_data, line


# Read element connectivities from interface file
def ReadElementConnectivities(interface_file,line,ElemsNr):
    # reading element connectivities
    elem_connectivities = np.zeros(4*ElemsNr, dtype=int)
    i = 0
    while line:
	if line == "\n":
            line = interface_file.readline()
        else:
            elem_connectivities[i*4+0] = int(line.split()[0])
            elem_connectivities[i*4+1] = int(line.split()[1])
            elem_connectivities[i*4+2] = int(line.split()[2])
            elem_connectivities[i*4+3] = int(line.split()[3])
            i = i+1
            line = interface_file.readline()

    return elem_connectivities


# Calculate cell pressure averaging nodal pressures
def CalculateCellPressure(nodal_pressures, node_ids):
    cell_pressure = 0.0

    # Interpolating nodal pressures
    for node in range(4):
        cell_pressure += 0.25 * nodal_pressures[node_ids[node]]

    return cell_pressure


# Calculate cell area



# Looks for a file matching the given pattern within the given path
def FindInitialMeshFilename(mesh_path, pattern):
    files_list = glob.glob(mesh_path + "*")
    for file in files_list:
        if file.startswith('%s' % mesh_path + '%s' % pattern):
            return file
    raise Exception('File: "{}" not found'.format(mesh_path + pattern))


# Looks for a file matching the given pattern within the given path
def FindMeshFilename(path, pattern, step):
    files_list = glob.glob(path + "*")
    if echo_level > 0:
        print('files_list = ', files_list)
    for file in files_list:
        if file.startswith('%s' % path + '%s' % pattern + '%s' % step):
            print(file)
            return file
    raise Exception('File: "{}" not found'.format(path + pattern + str(step)))


# Looks for a file matching the given pattern within the given path
def FindFilename(path, pattern, step):
    files_list = glob.glob(path + "*")
    #if echo_level > 0:
       ###08/07/2020      ###  print('files_list = ', files_list)
    for file in files_list:
        if file.startswith('%s' % path + '%s' % pattern + '%s' % step) and 'domain' not in file:
            print(file)
            return file
    raise Exception('File: "{}" not found'.format(path + pattern + str(step)))


# Calculate distance vector between two nodes
def CalculateDistanceVector(X,Y,Z,start,end):
    # Declaring and initializing distance vector
    distance_vector = np.zeros(3)

    # Computing distance
    distance_vector[0] = X[end] - X[start]
    distance_vector[1] = Y[end] - Y[start]
    distance_vector[2] = Z[end] - Z[start]

    return distance_vector



