from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# including kratos path
from KratosMultiphysics import *  # we import the KRATOS
from KratosMultiphysics.NurbsTestcaseApplication import * \
    # and now our application. note that we can import as many as we need to
    # solve our specific problem
from KratosMultiphysics.FluidDynamicsApplication import * \
    # and now our application. note that we can import as many as we need to
    # solve our specific problem
from nurbsdefinition import *
import time
# from line 13 to 35 a test fluid mesh can be created to test the distance computation of
# the NURBS-geometry. If commented out, the ModelPart will be created and stored in the
# working direction of this file.

import cube_mesher

xmin = 0
xmax = 50
ymin = 0
ymax = 50
zmin = 0
zmax = 50
nx = 10
ny = 10
nz = 10


box = cube_mesher.box_data(xmin, ymin, zmin, xmax, ymax, zmax, nx, ny, nz)

elemtype = "FractionalStepDiscontinuous3D"

with open("FluidMesh.mdpa", "w") as mdpa:
    cube_mesher.write_header(mdpa)
    cube_mesher.generate_nodes(mdpa, box)
    cube_mesher.generate_elements(mdpa, box, elemtype)
    cube_mesher.generate_mesh_groups(mdpa, box)

print("before read-in Fluid-Mesh")
print(time.clock())
# Read in the ModelPart "FluidMesh.mdpa" which has to be in the same
# working direction as this file
FluidMesh = ModelPart("FluidMesh")
mdpa_file = "FluidMesh"
model_part_io = ModelPartIO(mdpa_file)
model_part_io.ReadModelPart(FluidMesh)
print("after read-in Fluid-Mesh")
print(time.clock())

# crucial for solver (NURBS can so far only be computed in 2D)
domain_size = 2

NurbsModelPart = ModelPart("NurbsModelPart")
OutModelPart = ModelPart("OutModelPart")
                         # Only needed for visualization of
                         # Pure-Difusion-Problem

print("before read-in NurbsModelPart")
print(time.clock())
# Control Points are read in by mdpa-file
mdpa_file = "HalfSphere"
model_part_io = ModelPartIO(mdpa_file)
model_part_io.ReadModelPart(NurbsModelPart)
print("after read-in NurbsModelPart")
print(time.clock())

NurbsModeler = NurbsModeler()


print("before create half sphere bottom")
print(time.clock())

# Half Sphere bottom
polynomial_degree_xi = 2
polynomial_degree_eta = 2
knots_xi = [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0]
knots_eta = [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0]
number_of_cp_xi = 6
number_of_cp_eta = 6

weights = [
    1.0, 0.92677669529663687, 0.85355339059327373, 0.85355339059327395, 0.92677669529663698, 1.0,
    0.92677669529663687, 0.87185921676911449, 0.81694173824159222, 0.81694173824159222, 0.8718592167691146, 0.92677669529663687,
    0.85355339059327373, 0.81694173824159222, 0.78033008588991071, 0.78033008588991071, 0.81694173824159222, 0.85355339059327373,
    0.85355339059327395, 0.81694173824159222, 0.78033008588991071, 0.78033008588991071, 0.81694173824159222, 0.85355339059327395,
    0.92677669529663698, 0.8718592167691146, 0.81694173824159222, 0.81694173824159222, 0.8718592167691146, 0.92677669529663698,
    1.0, 0.92677669529663687, 0.85355339059327373, 0.85355339059327373, 0.92677669529663687, 1.0]

# connectivities state the first and last control point of a NURBS-Patch
connectivities_start = 1
connectivities_end = 36

# Using a python interface
NurbsSurface = NurbsDefinition(connectivities_start,
                               connectivities_end,
                               polynomial_degree_xi,
                               polynomial_degree_eta,
                               knots_xi,
                               knots_eta,
                               number_of_cp_xi,
                               number_of_cp_eta,
                               weights)

# Will create the elements of the patch
NurbsModeler.ReadModelPart(NurbsModelPart,
                           NurbsSurface.connectivities_start,
                           NurbsSurface.connectivities_end,
                           NurbsSurface.polynomial_degree_xi,
                           NurbsSurface.polynomial_degree_eta,
                           NurbsSurface.knots_xi,
                           NurbsSurface.knots_eta,
                           NurbsSurface.number_of_cp_xi,
                           NurbsSurface.number_of_cp_eta,
                           NurbsSurface.weights)


print("after create half sphere bottom")
print(time.clock())

#Half Sphere_ Quarter x+


print("before create half sphere Quarter_x+")
print(time.clock())
polynomial_degree_xi = 2
polynomial_degree_eta = 2
knots_xi = [0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0]
knots_eta = [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0]
number_of_cp_xi = 7
number_of_cp_eta = 6


weights = [1,		0.8535533906,	0.8535533906,	1,		0.8535533906,	0.8535533906,	1,
           0.9267766953,	0.7910533906,	0.7910533906,	0.9267766953,	0.7910533906,	0.7910533906,	0.9267766953,
           0.8535533906,	0.7285533906,	0.7285533906,	0.8535533906,	0.7285533906,	0.7285533906,	0.8535533906,
           0.8535533906,	0.7285533906,	0.7285533906,	0.8535533906,	0.7285533906,	0.7285533906,	0.8535533906,
           0.9267766953,	0.7910533906,	0.7910533906,	0.9267766953,	0.7910533906,	0.7910533906,	0.9267766953,
           1,		0.8535533906,	0.8535533906,	1,		0.8535533906,	0.8535533906,	1]


# connectivities state the first and last control point of a NURBS-Patch
connectivities_start = 37
connectivities_end = 78

# Using a python interface
NurbsSurface = NurbsDefinition(connectivities_start,
                               connectivities_end,
                               polynomial_degree_xi,
                               polynomial_degree_eta,
                               knots_xi,
                               knots_eta,
                               number_of_cp_xi,
                               number_of_cp_eta,
                               weights)

# Will create the elements of the patch

NurbsModeler.ReadModelPart(NurbsModelPart,
                           NurbsSurface.connectivities_start,
                           NurbsSurface.connectivities_end,
                           NurbsSurface.polynomial_degree_xi,
                           NurbsSurface.polynomial_degree_eta,
                           NurbsSurface.knots_xi,
                           NurbsSurface.knots_eta,
                           NurbsSurface.number_of_cp_xi,
                           NurbsSurface.number_of_cp_eta,
                           NurbsSurface.weights)


print("after create half sphere Quarter_x+")
print(time.clock())

# Half Sphere_ Quarter x-


print("before create half sphere Quarter_x-")
print(time.clock())
polynomial_degree_xi = 2
polynomial_degree_eta = 2
knots_xi = [0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0]
knots_eta = [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0]
number_of_cp_xi = 7
number_of_cp_eta = 6


weights = [	1,		0.8535533906,	0.8535533906,	1.0,		0.8535533906,	0.8535533906,
            1,		0.9267766953,	0.7910533906,	0.7910533906,	0.9267766953,	0.7910533906,
            0.7910533906,	0.9267766953,	0.8535533906,	0.7285533906,	0.7285533906,	0.8535533906,
            0.7285533906,	0.7285533906,	0.8535533906,	0.8535533906,	0.7285533906,	0.7285533906,
            0.8535533906,	0.7285533906,	0.7285533906,	0.8535533906,	0.9267766953,	0.7910533906,
            0.7910533906,	0.9267766953,	0.7910533906,	0.7910533906,	0.9267766953,	1,
            0.8535533906,	0.8535533906,	1,		0.8535533906,	0.8535533906,	1]


# connectivities state the first and last control point of a NURBS-Patch
connectivities_start = 79
connectivities_end = 120

# Using a python interface
NurbsSurface = NurbsDefinition(connectivities_start,
                               connectivities_end,
                               polynomial_degree_xi,
                               polynomial_degree_eta,
                               knots_xi,
                               knots_eta,
                               number_of_cp_xi,
                               number_of_cp_eta,
                               weights)

# Will create the elements of the patch
NurbsModeler.ReadModelPart(NurbsModelPart,
                           NurbsSurface.connectivities_start,
                           NurbsSurface.connectivities_end,
                           NurbsSurface.polynomial_degree_xi,
                           NurbsSurface.polynomial_degree_eta,
                           NurbsSurface.knots_xi,
                           NurbsSurface.knots_eta,
                           NurbsSurface.number_of_cp_xi,
                           NurbsSurface.number_of_cp_eta,
                           NurbsSurface.weights)
print("after create half sphere Quarter_x-")
print(time.clock())

# Only needed for Temperature Simulation
# import pure_diffusion_solver
# pure_diffusion_solver.AddVariables(NurbsModelPart)
# pure_diffusion_solver.AddVariables(OutModelPart)
# OutModelPart.AddNodalSolutionStepVariable(NURBS_ID)
# OutModelPart.AddNodalSolutionStepVariable(NURBS_COORDINATES)

# mdpa_file = "Temperature_Example"
# model_part_io = ModelPartIO(mdpa_file)
# model_part_io.ReadModelPart(OutModelPart)


# Distance and Radius are needed for the Closest-Point-Search
FluidMesh.AddNodalSolutionStepVariable(DISTANCE)
FluidMesh.AddNodalSolutionStepVariable(RADIUS)


print("before calling ClosestPoint")
print(time.clock())
# running the Closes-Point-Search
Max_search_radius = 20.0
NurbsModeler.ClosestPoint(FluidMesh, NurbsModelPart, Max_search_radius)
print("after calling ClosestPoint")
print(time.clock())


# Boundary Conditions for Temperature Problem

# NurbsModelPart.SetBufferSize(1)
# pure_diffusion_solver.AddDofs(NurbsModelPart)
# 2D Temperature Problem
# for node in NurbsModelPart.Nodes:
#    if node.X == 1:
#        node.Fix(TEMPERATURE)
#        node.SetSolutionStepValue(TEMPERATURE,0,10.0)
#    elif node.Y == 2:
#        node.Fix(TEMPERATURE)
#        node.SetSolutionStepValue(TEMPERATURE,0,0)


# Unit Square 2 Elements
# for node in NurbsModelPart.Nodes:
#        if node.X == 0:
#            node.Fix(TEMPERATURE)
#            node.SetSolutionStepValue(TEMPERATURE,0,5.0)
#        elif node.X == 1:
#            node.Fix(TEMPERATURE)
#            node.SetSolutionStepValue(TEMPERATURE,0,0)
#        elif node.X == 0.5:
#            node.SetSolutionStepValue(TEMPERATURE,0,2.5)
#    node.SetSolutionStepValue(TEMPERATURE,0,5.0*(1.0-node.X))

# solver = pure_diffusion_solver.StaticPoissonSolver(NurbsModelPart,domain_size)
# solver.Initialize()
# solver.SetEchoLevel(3)
# solver.Solve()


# now we proceed to use the GID interface (both to import the infomation
# inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostAscii  # PostBinary for getting readable info  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("control_points", gid_mode,
               multifile, deformed_mesh_flag, write_conditions)

# Improvised Post-Processing with TriangleMesh
# NurbsModeler.InterpolateDesignVariables(NurbsModelPart,OutModelPart)

# we create a mesh for the postprocess
# mesh_name = 0.0
# gid_io.InitializeMesh( mesh_name );
# gid_io.WriteNodeMesh((NurbsModelPart).GetMesh());
# gid_io.FinalizeMesh()

# gid_io.InitializeResults(mesh_name,(NurbsModelPart).GetMesh())
# gid_io.WriteNodalResults(TEMPERATURE,NurbsModelPart.Nodes,0,0)
# gid_io.FinalizeResults()

# gid_io_out = GidIO("triange_visualization",gid_mode,multifile,deformed_mesh_flag,write_conditions)
# mesh_name = 0.0
# gid_io_out.InitializeMesh( mesh_name );
# gid_io_out.WriteMesh((OutModelPart).GetMesh());
# gid_io_out.FinalizeMesh()

# gid_io_out.InitializeResults(mesh_name,(OutModelPart).GetMesh())
# gid_io_out.WriteNodalResults(NURBS_COORDINATES,OutModelPart.Nodes,0,0)
# gid_io_out.WriteNodalResults(TEMPERATURE,OutModelPart.Nodes,0,0)
# gid_io_out.FinalizeResults()


# Write Fluidmesh for visualizing the Solutions of the Closes-Point-Search
mesh_name = 0.0
gid_io_out = GidIO("Closest_Point_Search", gid_mode,
                   multifile, deformed_mesh_flag, write_conditions)
gid_io_out.InitializeMesh(mesh_name)
gid_io_out.WriteMesh((FluidMesh).GetMesh())
gid_io_out.FinalizeMesh()

gid_io_out.InitializeResults(mesh_name, (FluidMesh).GetMesh())
gid_io_out.WriteNodalResults(DISTANCE, FluidMesh.Nodes, 0, 0)
gid_io_out.WriteNodalResults(RADIUS, FluidMesh.Nodes, 0, 0)
gid_io_out.FinalizeResults()
