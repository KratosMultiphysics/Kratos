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


# crucial for solver (NURBS can so far only be computed in 2D)
domain_size = 2

NurbsModelPart = ModelPart("NurbsModelPart")
OutModelPart = ModelPart("OutModelPart")

import pure_diffusion_solver
pure_diffusion_solver.AddVariables(NurbsModelPart)
pure_diffusion_solver.AddVariables(OutModelPart)

# Control Points are read in by mdpa-file
mdpa_file = "2DTemperatureProblemKratosDummies"
print(mdpa_file)

NurbsModelPart.AddNodalSolutionStepVariable(TEMPERATURE)

OutModelPart.AddNodalSolutionStepVariable(NURBS_ID)
OutModelPart.AddNodalSolutionStepVariable(NURBS_COORDINATES)

model_part_io = ModelPartIO(mdpa_file)
model_part_io.ReadModelPart(NurbsModelPart)
NurbsModelPart.SetBufferSize(2)

NurbsModeler = NurbsModeler()


# Temperature Example


polynomial_degree_xi = 3
polynomial_degree_eta = 1
knots_xi = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
knots_eta = [0.0, 0.0, 1.0, 1.0]
number_of_cp_xi = 5
number_of_cp_eta = 2


weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# connectivities state the first and last control point of a NURBS-Patch
connectivities_start = 1
connectivities_end = 10

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
                           domain_size,
                           NurbsSurface.connectivities_start,
                           NurbsSurface.connectivities_end,
                           NurbsSurface.polynomial_degree_xi,
                           NurbsSurface.polynomial_degree_eta,
                           NurbsSurface.knots_xi,
                           NurbsSurface.knots_eta,
                           NurbsSurface.number_of_cp_xi,
                           NurbsSurface.number_of_cp_eta,
                           NurbsSurface.weights)

# Only needed for Temperature Simulation



mdpa_file = "Temperature_Example"
model_part_io = ModelPartIO(mdpa_file)
model_part_io.ReadModelPart(OutModelPart)


# Boundary Conditions for Temperature Problem
NurbsModelPart.SetBufferSize(2)

pure_diffusion_solver.AddDofs(NurbsModelPart)
# 2D Temperature Problem
for node in NurbsModelPart.Nodes:
    if node.X == 1:
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE, 0, 10.0)
    elif node.Y == 2:
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE, 0, 0)


solver = pure_diffusion_solver.StaticPoissonSolver(NurbsModelPart, domain_size)
solver.Initialize()
solver.SetEchoLevel(3)
solver.Solve()


# now we proceed to use the GID interface (both to import the infomation
# inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostAscii  # PostBinary for getting readable info  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("control_points", gid_mode,
               multifile, deformed_mesh_flag, write_conditions)

# Improvised Post-Processing with TriangleMesh
NurbsModeler.InterpolateDesignVariables(NurbsModelPart, OutModelPart)

# we create a mesh for the postprocess
mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteNodeMesh((NurbsModelPart).GetMesh())
gid_io.FinalizeMesh()

gid_io.InitializeResults(mesh_name, (NurbsModelPart).GetMesh())
gid_io.WriteNodalResults(TEMPERATURE, NurbsModelPart.Nodes, 0, 0)
gid_io.FinalizeResults()

gid_io_out = GidIO("triange_visualization", gid_mode,
                   multifile, deformed_mesh_flag, write_conditions)
mesh_name = 0.0
gid_io_out.InitializeMesh(mesh_name)
gid_io_out.WriteMesh((OutModelPart).GetMesh())
gid_io_out.FinalizeMesh()

gid_io_out.InitializeResults(mesh_name, (OutModelPart).GetMesh())
gid_io_out.WriteNodalResults(NURBS_COORDINATES, OutModelPart.Nodes, 0, 0)
gid_io_out.WriteNodalResults(TEMPERATURE, OutModelPart.Nodes, 0, 0)
gid_io_out.FinalizeResults()


# Write Fluidmesh for visualizing the Solutions of the Closes-Point-Search
# mesh_name = 0.0
# gid_io_out = GidIO("Closest_Point_Search",gid_mode,multifile,deformed_mesh_flag,write_conditions)
# gid_io_out.InitializeMesh( mesh_name );
# gid_io_out.WriteMesh((FluidMesh).GetMesh());
# gid_io_out.FinalizeMesh()

# gid_io_out.InitializeResults(mesh_name,(FluidMesh).GetMesh())
# gid_io_out.WriteNodalResults(DISTANCE,FluidMesh.Nodes,0,0)
# gid_io_out.WriteNodalResults(RADIUS,FluidMesh.Nodes,0,0)
# gid_io_out.FinalizeResults()
