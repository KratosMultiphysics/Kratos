# -*- coding: utf-8 -*-
import re, glob, subprocess, time, os
import numpy as np
import tau_python
from tau_python import tau_msg
# import tau_python.tau_msg as tau_msg
import PyPara, PySurfDeflect
from scipy.io import netcdf

echo_level = 1

# GetFluidMesh is called only once at the beginning, after the first fluid solve
def GetFluidMesh(working_path, step, para_path_mod):
    # Find the interface file name
    interface_file_name = FindInterfaceFile(working_path, step)

    # Read mesh from interface file
    NodesNr, ElemsNr, X, Y, Z, P, elem_connectivities = ReadTauOutput(interface_file_name, 20)

    # Transform nodal coordinates to numpy array
    nodal_coords = ReadNodalCoordinates(NodesNr, X, Y, Z)

    # Save element types in a numpy array
    element_types = ReadElementTypes(ElemsNr)

    # In vtk format element connectivities start from 0, not from 1
    elem_connectivities -= 1

    return nodal_coords, elem_connectivities, element_types


def ComputeFluidForces(working_path, step):
    # Find the interface file name
    interface_file_name = FindInterfaceFile(working_path, step)

    # Read mesh and pressure from interface file
    NodesNr, ElemsNr, X, Y, Z, P, elem_connectivities = ReadTauOutput(interface_file_name, 20)

    # Calculate the cell pressure averaging the nodal pressure
    cell_pressure = CalculateCellPressure(ElemsNr, P, X, elem_connectivities)

    # Calculate cells' normals and areas
    cells_normals = CalculateCellNormals(ElemsNr, elem_connectivities, X, Y, Z)
    cells_areas = CalculateCellAreas(ElemsNr, elem_connectivities, X, Y, Z)

    # calculating the force vector
    fluid_forces = CalculateNodalFluidForces(
        ElemsNr, elem_connectivities, NodesNr, cell_pressure, cells_areas, cells_normals)

    return fluid_forces

def ExecuteBeforeMeshDeformation(dispTau, working_path, step, para_path_mod, start_step):
    global dispTauOld
    print "deformationstart"
    interface_file_name = FindInterfaceFile(working_path, step)

    NodesNr,ElemsNr,X,Y,Z,P,elemTable=ReadTauOutput(interface_file_name, 20)

    nodes = ReadNodalCoordinates(NodesNr, X, Y, Z)

    if(step==start_step):
        dispTauOld=np.zeros(3*NodesNr)
        dispTau_transpose = np.transpose(dispTau)
        print 'dispTau =', dispTau_transpose
    print 'dispTauOld = ', dispTauOld

    [ids,coordinates,globalID,coords]=meshDeformation(NodesNr,nodes,dispTau,dispTauOld,para_path_mod)
    PySurfDeflect.write_test_surface_file('deformation_file',coords[:,0:2],coords[:,3:5])
    print "afterPySurfDeflect"

    for i in xrange(0,3*NodesNr):
        dispTauOld[i]=dispTau[i]
    print "afterDeformation"


def FindOutputFile(working_path, step):
    outputs_path = working_path + "Outputs/"
    CheckIfPathExists(outputs_path)
    ouput_file_pattern = "airfoilSol.pval.unsteady_i="
    return FindFileName(outputs_path, ouput_file_pattern, step + 1)


def FindMeshFile(working_path, step, start_step):
    mesh_path = working_path + "Mesh/"
    CheckIfPathExists(mesh_path)
    if step == start_step:
        pattern = 'airfoil_Structured_scaliert.grid'
        return FindInitialMeshFileName(mesh_path, pattern)
    else:
        pattern = 'airfoil_Structured_scaliert.grid.def.'
        return FindFileName(mesh_path, pattern, step-start_step)


def ConvertOutputToDat(working_path, tau_path, step, para_path_mod, start_step):
    PrintBlockHeader("Start Writting Solution Data at time %s" % (str(time)))
    subprocess.call('rm ' + working_path + '/Tautoplt.cntl', shell=True)
    tautoplt_file_name = WriteTautoplt(working_path, step, para_path_mod, start_step)
    command = tau_path + 'tau2plt ' + tautoplt_file_name
    subprocess.call(command, shell=True)
    PrintBlockHeader("Stop Writting Solution Data at time %s" % (str(time)))


def FindInterfaceFile(working_path, step):
    output_file_name = FindOutputFile(working_path, step)
    interface_file_name = output_file_name.replace('pval', 'surface.pval')
    interface_file_name = interface_file_name.replace('+', '') + '.dat'
    CheckIfPathExists(interface_file_name)
    return interface_file_name


def CheckIfPathExists(path):
    if not os.path.exists(path):
        raise Exception('Path: "{}" not found'.format(path))


def FindInitialMeshFileName(path, name):
    files_list = glob.glob(path + "*")
    for file in files_list:
        if file.startswith('%s' % path + '%s' % name):
            return file
    raise Exception('File: "{}" not found'.format(path + name))


def FindFileName(path, name, step):
    files_list = glob.glob(path + "*")
    if echo_level > 0:
        print 'files_list =', files_list
    for file in files_list:
        if file.startswith('%s' % path + '%s' % name + '%s' % step):
            return file
    raise Exception('File: "{}" not found'.format(path + name + str(step)))


def WriteTautoplt(working_path, step, para_path_mod, start_step):
    mesh_file_name = FindMeshFile(working_path, step, start_step)
    parameter_file_name = working_path + para_path_mod
    output_file_name = FindOutputFile(working_path, step)
    tautoplt_file_name = working_path + 'Tautoplt.cntl'
    tautoplt_file_writing = open(tautoplt_file_name, 'w')
    tautoplt_file_reading = open(working_path + 'Tautoplt_initial.cntl', 'r+')
    line = tautoplt_file_reading.readline()
    while line:
        if 'Primary grid filename:' in line:
            line = 'Primary grid filename:' + mesh_file_name + ' \n'
            tautoplt_file_writing.write(line)
            line = tautoplt_file_reading.readline()
        if 'Boundary mapping filename:' in line:
            line = 'Boundary mapping filename:' + parameter_file_name + ' \n'
            tautoplt_file_writing.write(line)
            line = tautoplt_file_reading.readline()
        if 'Restart-data prefix:' in line:
            line = 'Restart-data prefix:' + output_file_name + ' \n'
            tautoplt_file_writing.write(line)
            line = tautoplt_file_reading.readline()
        else:
            line = tautoplt_file_reading.readline()
            tautoplt_file_writing.write(line)
    tautoplt_file_writing.close()
    tautoplt_file_reading.close()
    return tautoplt_file_name



def PrintBlockHeader(header):
    tau_python.tau_msg("\n" + 50 * "*" + "\n" + "* %s\n" %header + 50*"*" + "\n")

# Read mesh and data from tau output file
def ReadTauOutput(interface_file_name, velocity):
    position_info, mesh_info, nodal_data, elem_connectivities = ReadInterfaceFile(
        interface_file_name)

    # Read mesh info
    NodesNr = mesh_info[0]
    ElemsNr = mesh_info[1]

    X, Y, Z = SaveCoordinatesList(nodal_data, position_info, NodesNr)
    P = SavePressure(nodal_data, position_info, NodesNr, velocity)

    return NodesNr, ElemsNr, X, Y, Z, P, elem_connectivities

def ReadInterfaceFile(interface_file_name):
    with open(interface_file_name, 'r') as interface_file:
        line = interface_file.readline()
        position_info, mesh_info, line = ReadHeader(interface_file, line)
        nodal_data, line = ReadNodalData(interface_file, line)
        elem_connectivities = ReadElementConnectivities(interface_file, line, ElemsNr=mesh_info[1])
    return position_info, mesh_info, nodal_data, elem_connectivities

def SaveCoordinatesList(nodal_data, position_info, NodesNr):
    # Read coordinates positions and save them in a list
    x_position = position_info.index('"x"') - 2
    X = nodal_data[(x_position)*NodesNr:(x_position+1)*NodesNr]
    y_position = position_info.index('"y"') - 2
    Y = nodal_data[(y_position)*NodesNr:(y_position+1)*NodesNr]
    z_position = position_info.index('"z"') - 2
    Z = nodal_data[(z_position)*NodesNr:(z_position+1)*NodesNr]
    return X, Y, Z

def SavePressure(nodal_data, position_info, NodesNr, velocity):
    # Read pressure coefficient
    cp_position = position_info.index('"cp"') - 2
    CP = nodal_data[(cp_position)*NodesNr:(cp_position+1)*NodesNr]
    # Compute pressure from cp
    P = [x*velocity*velocity*0.5*1.2 for x in CP]
    P = np.array(P)
    return P

def ReadNodalData(interface_file,line):
    # reading nodal data
    nodal_data = []
    while 'E+' in line or 'E-' in line:
        for elem in line.split():
            nodal_data.append(float(elem))
        line = interface_file.readline()

    return nodal_data, line

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


def ReadNodalCoordinates(NodesNr, X, Y, Z):
    # array to store the coordinates of the nodes: x1,y1,z1,x2,y2,z2,...
    nodal_coords = np.zeros(NodesNr*3)

    for i in xrange(0, NodesNr):
        nodal_coords[3*i+0] = X[i]
        nodal_coords[3*i+1] = Y[i]
        nodal_coords[3*i+2] = Z[i]

    return nodal_coords


def ReadElementTypes(ElemsNr):
    # array to store the element types
    element_types = np.zeros(ElemsNr, dtype=int)

    for i in xrange(0, ElemsNr):
        element_types[i] = 9

    return element_types

def ReadElemConnectivities(ElemsNr, elemTable):
    # array to store the element connectivities
    elem_connectivities = np.zeros(4*ElemsNr, dtype=int)
    element_types = np.zeros(ElemsNr, dtype=int)

    for i in xrange(0, ElemsNr):
        elem_connectivities[i*4+0] = elemTable[i, 0]
        elem_connectivities[i*4+1] = elemTable[i, 1]
        elem_connectivities[i*4+2] = elemTable[i, 2]
        elem_connectivities[i*4+3] = elemTable[i, 3]
        element_types[i] = 9

    return elem_connectivities, element_types


# Calculate the cell pressure averaging the nodal pressure
def CalculateCellPressure(ElemsNr, P, X, elem_connectivities):
    cell_pressure = np.zeros(ElemsNr)  # cp for interface elements

    for i in range(ElemsNr):
        for k in range(4):
            cell_pressure[i] += 0.25 * P[elem_connectivities[i*4+k]-1]

    if echo_level > 0:
        with open('xp', 'w') as f:
            for i in range(ElemsNr):
                x = 0.0
                for k in range(4):
                    x += 0.25 * X[elem_connectivities[i*4+k]-1]
                f.write('%d\t%f\t%f\n' % (i, x, cell_pressure[i]))

    return cell_pressure


def CalculateCellNormals(ElemsNr, elem_connectivities, X, Y, Z):
    cells_normals = np.zeros([ElemsNr, 3])

    # Loop over all cells
    for i in xrange(0, ElemsNr):
        node_ids = GetCellNodeIds(elem_connectivities, i)
        # Calculate cell diagonals
        cell_diagonal_02 = CalculateDistanceVector(X,Y,Z,node_ids[0],node_ids[2])
        cell_diagonal_13 = CalculateDistanceVector(X,Y,Z,node_ids[1],node_ids[3])

        # Calculate cell normal
        cell_normal = np.cross(cell_diagonal_02, cell_diagonal_13)
        magnitude = np.linalg.norm(cell_normal)
        unit_cell_normal = cell_normal/magnitude
        cells_normals[i, 0] = unit_cell_normal[0]
        cells_normals[i, 1] = unit_cell_normal[1]
        cells_normals[i, 2] = unit_cell_normal[2]

    return cells_normals


def CalculateCellAreas(ElemsNr, elem_connectivities, X, Y, Z):
    cells_areas = np.zeros(ElemsNr)

    # Loop over all cells
    for i in xrange(0, ElemsNr):
        node_ids = GetCellNodeIds(elem_connectivities, i)

        # Calculate cell sides
        cell_side_01 = CalculateDistanceVector(X, Y, Z, node_ids[0], node_ids[1])
        cell_side_03 = CalculateDistanceVector(X, Y, Z, node_ids[0], node_ids[3])

        # Calculate cell area
        cell_area = np.cross(cell_side_01, cell_side_03)

        # hold only for 2d case, MUST be modified for 3d interfaces
        cells_areas[i] = np.linalg.norm(cell_area)

    return cells_areas

def GetCellNodeIds(elem_connectivities, elem_id):
    node_ids = np.zeros(4, dtype=int)
    node_ids[0] = elem_connectivities[elem_id*4+0]-1
    node_ids[1] = elem_connectivities[elem_id*4+1]-1
    node_ids[2] = elem_connectivities[elem_id*4+2]-1
    node_ids[3] = elem_connectivities[elem_id*4+3]-1
    return node_ids

def CalculateDistanceVector(X,Y,Z,start,end):
    distance_vector = np.zeros(3)
    distance_vector[0] = X[end] - X[start]
    distance_vector[1] = Y[end] - Y[start]
    distance_vector[2] = Z[end] - Z[start]
    return distance_vector

def CalculateCellForce(pressure,area,normal,cell_id):
    cell_force = np.zeros(3)
    cell_force[0] = pressure * area * normal[cell_id, 0]
    cell_force[1] = pressure * area * normal[cell_id, 1]
    cell_force[2] = pressure * area * normal[cell_id, 2]
    return cell_force


# Calculate the Vector Force
def CalculateNodalFluidForces(ElemsNr,elem_connectivities,NodesNr,cell_pressure,area,normal):
    forcesTauNP = np.zeros(NodesNr*3)
    for i in xrange(0, ElemsNr):
        cell_force = CalculateCellForce(cell_pressure[i],area[i],normal,i)
        for k in range(4):
            for j in range(3):
                forcesTauNP[3*(elem_connectivities[i*4+k]-1)+j] += 0.25 * cell_force[j]

    return forcesTauNP

# Execute the Mesh deformation of TAU
def meshDeformation(NodesNr,nodes,dispTau,dispTauOld, para_path_mod):

    disp=np.zeros([NodesNr,3])#NodesNr
    for i in xrange(0,NodesNr):#NodesNr
        disp[i,0]=1*(dispTau[3*i+0]-dispTauOld[3*i+0])
        disp[i,1]=1*(dispTau[3*i+1]-dispTauOld[3*i+1])
        disp[i,2]=1*(dispTau[3*i+2]-dispTauOld[3*i+2])
    Para = PyPara.Parafile(para_path_mod)
    ids, coordinates = PySurfDeflect.read_tau_grid(Para)
    coords=np.zeros([NodesNr,6])#NodesNr

    for i in xrange(0,NodesNr):
        coords[i,0]=coordinates[0,i]
        coords[i,1]=coordinates[1,i]
        coords[i,2]=coordinates[2,i]

    globalID = np.zeros(NodesNr)
    for i in xrange(0,NodesNr):
        xi = coords[i,0]
        yi = coords[i,1]
        zi = coords[i,2]
        for k in xrange(0,NodesNr):
            xk = nodes[3*k+0]
            yk = nodes[3*k+1]
            zk = nodes[3*k+2]
            dist2 = (xi-xk)*(xi-xk) + (yi-yk)*(yi-yk) + (zi-zk)*(zi-zk)

            if dist2 < 0.00001:
                #K=k
                globalID[i] = k
            globalID = globalID.astype(int)
        #print "%d found %d" % (i,K)
        coords[i, 3] = disp[globalID[i], 0]
        coords[i, 4] = disp[globalID[i], 1]
        coords[i, 5] = disp[globalID[i], 2]

    fname_new = 'interface_deformfile.nc'
    ncf = netcdf.netcdf_file(fname_new, 'w')
    # define dimensions
    nops = 'no_of_points'
    number_of_points = len(ids[:])
    ncf.createDimension(nops, number_of_points)
    # define variables
    gid = ncf.createVariable('global_id', 'i', (nops,))
    ncx = ncf.createVariable('x', 'd', (nops,))
    ncy = ncf.createVariable('y', 'd', (nops,))
    ncz = ncf.createVariable('z', 'd', (nops,))
    ncdx = ncf.createVariable('dx', 'd', (nops,))
    ncdy = ncf.createVariable('dy', 'd', (nops,))
    ncdz = ncf.createVariable('dz', 'd', (nops,))
    # write data
    gid[:] = ids
    ncx[:] = coords[:,0]
    ncy[:] = coords[:,1]
    ncz[:] = coords[:,2]
    ncdx[:] = coords[:,3]
    ncdy[:] = coords[:,4]
    ncdz[:] = coords[:,5]
    ncf.close()

    return ids, coordinates, globalID, coords
