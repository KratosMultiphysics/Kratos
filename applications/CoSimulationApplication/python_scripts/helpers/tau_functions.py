# -*- coding: utf-8 -*-
import re, glob, subprocess, time, os, warnings
import numpy as np
import tau_python
from tau_python import tau_msg
from tau_python import tau_solver_unsteady_get_physical_time
import PyPara, PySurfDeflect
from scipy.io import netcdf

echo_level = 1


# Convert tau output to dat file using tau2plt
def ConvertOutputToDat(working_path, tau_path, step, para_path_mod, start_step):
    PrintBlockHeader("Start Writting Solution Data at time %s" % (str(time)))

    # Write Tautoplt.cntl file
    tautoplt_filename = WriteTautoplt(working_path, step, para_path_mod, start_step)

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

    # Read tau's parameter file
    Para = PyPara.Parafile(para_path_mod)

    # Read the interface fluid grid
    ids, coordinates = PySurfDeflect.read_tau_grid(Para)

    # Write membrane's displacments in a file
    WriteInterfaceDeformationFile(ids, coordinates, relative_displacements)


# Computes fluid forces at the nodes
def ComputeFluidForces(working_path, step):
    # Read mesh and pressure from interface file
    X, Y, Z, nodal_pressures, elem_connectivities = ReadTauOutput(working_path, step, 20)

    # calculating the force vector
    fluid_forces = CalculateNodalFluidForces(X, Y, Z, nodal_pressures, elem_connectivities)

    return fluid_forces


# GetFluidMesh is called only once at the beginning, after the first fluid solve
def GetFluidMesh(working_path, step, para_path_mod):
    # Read mesh from interface file
    X, Y, Z, P, elem_connectivities = ReadTauOutput(working_path, step, 20)

    # Transform nodal coordinates to numpy array
    nodal_coords = ReadNodalCoordinates(X, Y, Z)

    # Save element types in a numpy array
    element_types = ReadElementTypes(len(elem_connectivities)/4)

    # In vtk format element connectivities start from 0, not from 1
    elem_connectivities -= 1

    return nodal_coords, elem_connectivities, element_types


# Write Tautoplt.cntl file
def WriteTautoplt(working_path, step, para_path_mod, start_step):
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
            line = ModifyFilesIOLines(line, working_path, step, para_path_mod, start_step)
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
    number_of_nodes = len(total_displacements)/3
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


# Write membrane's displacments in a file
def WriteInterfaceDeformationFile(ids, coordinates, relative_displacements):
    # Open interface_deformfile
    ncf = netcdf.netcdf_file('interface_deformfile.nc', 'w')

    # define dimensions
    nops = 'no_of_points'
    number_of_points = len(ids[:])
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


# Read mesh and data from tau output file
def ReadTauOutput(working_path, step, velocity):
    # Find the interface file name
    interface_filename = FindInterfaceFilename(working_path, step)

    # Read interface file
    position_info, mesh_info, nodal_data, elem_connectivities = ReadInterfaceFile(
        interface_filename)

    # Read mesh info
    NodesNr = mesh_info[0]
    ElemsNr = mesh_info[1]

    X, Y, Z = SaveCoordinatesList(nodal_data, position_info, NodesNr)
    P = SavePressure(nodal_data, position_info, NodesNr, velocity)

    return X, Y, Z, P, elem_connectivities


# Calculate the fluid forces at the nodes
def CalculateNodalFluidForces(X, Y, Z, nodal_pressures, elem_connectivities):
    nodal_forces = np.zeros(3*len(X))
    # Loop over cells
    for cell in range(len(elem_connectivities)/4):
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
    # array to store the element types
    element_types = np.zeros(ElemsNr, dtype=int)

    for i in xrange(0, ElemsNr):
        element_types[i] = 9

    return element_types


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
def ModifyFilesIOLines(line, working_path, step, para_path_mod, start_step):
    if 'Primary grid filename:' in line:
        primary_grid_filename = FindPrimaryGridFilename(working_path, step, start_step)
        line = 'Primary grid filename:' + primary_grid_filename + ' \n'
    elif 'Boundary mapping filename:' in line:
        parameter_filename = working_path + para_path_mod
        line = 'Boundary mapping filename:' + parameter_filename + ' \n'
    elif 'Restart-data prefix:' in line:
        output_filename = FindOutputFilename(working_path, step)
        line = 'Restart-data prefix:' + output_filename + ' \n'

    return line


# Find the interface file name
def FindInterfaceFilename(working_path, step):
    output_filename = FindOutputFilename(working_path, step)
    interface_filename = output_filename.replace('pval', 'surface.pval')
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
    P = [x*velocity*velocity*0.5*1.2 for x in CP]
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


# Finds steps' corresponding primary grid filename
def FindPrimaryGridFilename(working_path, step, start_step):
    mesh_path = working_path + "Mesh/"
    CheckIfPathExists(mesh_path)
    if step == start_step:
        pattern = 'airfoil_Structured_scaliert.grid'
        return FindInitialMeshFilename(mesh_path, pattern)
    else:
        pattern = 'airfoil_Structured_scaliert.grid.def.'
        return FindFilename(mesh_path, pattern, step-start_step)


# Finds output filename
def FindOutputFilename(working_path, step):
    outputs_path = working_path + "Outputs/"
    CheckIfPathExists(outputs_path)
    ouput_file_pattern = "airfoilSol.pval.unsteady_i="
    return FindFilename(outputs_path, ouput_file_pattern, step + 1)


# Read header from interface file
def ReadHeader(interface_file,line):
    # reading the five lines of the header
    for _ in range(5):
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


# Looks for a file matching the given pattern within the given path
def FindInitialMeshFilename(mesh_path, pattern):
    files_list = glob.glob(mesh_path + "*")
    for file in files_list:
        if file.startswith('%s' % mesh_path + '%s' % pattern):
            return file
    raise Exception('File: "{}" not found'.format(mesh_path + pattern))


# Looks for a file matching the given pattern within the given path
def FindFilename(path, pattern, step):
    files_list = glob.glob(path + "*")
    if echo_level > 0:
        print 'files_list =', files_list
    for file in files_list:
        if file.startswith('%s' % path + '%s' % pattern + '%s' % step):
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

class MotionStringGenerator(object):
	"""Auxiliary class to generate TAU motion strings for a translatory x-motion
	and/or a pitching oscillation.
	"""

	def __init__(self, deltaT, pitchDeg, thetaDeg, thetaRate):
		self.deltaT       = deltaT
		self.pitchDeg = pitchDeg
		self.thetaDeg = thetaDeg
		self.thetaRate = thetaRate
		

	def GetMotionString(self,step):
		self.time = step*self.deltaT
		self.thetaInstant     = self.thetaDeg[step]
		self.pitchFreq    = self.thetaRate[step]

		p     = 0.
		q     = np.deg2rad(self.pitchFreq)
		r     = 0.
		phi   = 0.
		theta = np.deg2rad(self.pitchDeg) + np.deg2rad(self.thetaInstant)
		psi   = 0.
		u     = 0
		v     = 0.
		w     = 0.
		dx    = 0.
		dy    = 0.
		dz    = 0.
		motionString=" ".join(map(str, [p,q,r,phi,theta,psi,u,v,w,dx,dy,dz]))

		return motionString

	def __call__(self, step):
		return self.GetMotionString(step)
