from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# setting the domain size for the problem to be solved
domain_size = 2

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.ALEApplication import *

# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
structure_model_part = ModelPart("StructurePart")

#
# importing the solvers needed
import mesh_solver
import incompressible_fluid_solver
# import NonConformant_OneSideMap
import Conformant_OneSideMap
import structural_solver_dynamic
import FractionalStepCouplingAitken

incompressible_fluid_solver.AddVariables(fluid_model_part)
mesh_solver.AddVariables(fluid_model_part)
# NonConformant_OneSideMap.AddVariables(fluid_model_part,structure_model_part)
Conformant_OneSideMap.AddVariables(fluid_model_part, structure_model_part)
structural_solver_dynamic.AddVariables(structure_model_part)
FractionalStepCouplingAitken.AddVariables(fluid_model_part, structure_model_part)

# introducing input file name
input_file_name = "flag_test"

# reading the fluid part
gid_io = GidIO(input_file_name + str("_fluid"), GiDPostMode.GiD_PostBinary)
gid_io.ReadModelPart(fluid_model_part)
print(fluid_model_part)
print("fluid model read correctly")

# reading the structural part
data_io = DatafileIO(input_file_name + str("_structure"))
data_io.ReadModelPart(structure_model_part)
print(structure_model_part)
print("structural model read correctly")

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)
structure_model_part.SetBufferSize(3)

# adding dofs
incompressible_fluid_solver.AddDofs(fluid_model_part)
mesh_solver.AddDofs(fluid_model_part)
structural_solver_dynamic.AddDofs(structure_model_part)


# assigning the fluid properties
density = 1.18
viscosity = 1.82e-5 / density
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY, 0, viscosity)
    node.SetSolutionStepValue(DENSITY, 0, density)


# creating the solvers
# fluid solver
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part, domain_size)
fluid_solver.laplacian_form = 1
# standard laplacian form
fluid_solver.predictor_corrector = True
fluid_solver.max_press_its = 10
fluid_solver.Initialize()
print("fluid solver created")

# mesh solver
reform_dofs_at_each_step = False
mesh_solver = mesh_solver.MeshSolver(fluid_model_part, domain_size, reform_dofs_at_each_step)
pDiagPrecond = DiagonalPreconditioner()
mesh_solver.linear_solver = CGSolver(1e-3, 300, pDiagPrecond)
mesh_solver.time_order = 2
mesh_solver.Initialize()
mesh_solver.solver.SetEchoLevel(0)
print("mesh solver created")

# structure solver
structure_solver = structural_solver_dynamic.DynamicStructuralSolver(structure_model_part, domain_size)
structure_solver.echo_level = 0
# structure_solver.toll = 1e-12
# structure_solver.absolute_tol = 1e-15

# pILUPrecond = ILU0Preconditioner()
# structure_solver.structure_linear_solver = BICGSTABSolver(1e-9, 5000,pILUPrecond)
structure_solver.Initialize()
structure_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D())
print("structural solver created")

# mapper
# non conformant mapper
# mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(fluid_model_part,structure_model_part)
# conformant point to point
utilities = VariableUtils()
interface_fluid_nodes = (utilities).SelectNodeList(IS_INTERFACE, 1.0, fluid_model_part.Nodes)
interface_structure_nodes = (utilities).SelectNodeList(IS_INTERFACE, 1.0, structure_model_part.Nodes)
print("interface fluid nodes = ", len(interface_fluid_nodes))
print("interface structure nodes = ", len(interface_structure_nodes))
mapper = Conformant_OneSideMap.Conformant_OneSideMap(interface_fluid_nodes, interface_structure_nodes)

print("mapper created")

# creating the coupled solver
coupled_solver = FractionalStepCouplingAitken.FractionalStepCoupling(fluid_model_part, structure_model_part, structure_solver, mesh_solver, mapper, domain_size)
coupled_solver.Initialize()

print("coupled solver created")


def FindNode(node_list, x, y, z):
    for node in node_list:
        a = (node.X - x) ** 2
        a = a + (node.Y - y) * (node.Y - y)
        a = a + (node.Z - z) * (node.Z - z)
        if (a < 0.0000001):
            return node


def PrintDisp(time, filename, node):
    print("in print force", time)
    outstring = str(time) + " "
    outstring += str(node.GetSolutionStepValue(DISPLACEMENT_Y)) + "\n"
    filename.write(outstring)
    filename.flush()

tip = FindNode(structure_model_part.Nodes, 0.105, 0.06, 0.0)

scalarout = open("tip_disp.csv", 'w')
scalarout.write("time tip_disp \n")


# settings to be changed
nsteps = 5000
output_step = 5

Dt = 0.005

out = 0

# mesh to be printed
gid_io.WriteMesh((fluid_model_part).GetMesh(), domain_size, GiDPostMode.GiD_PostBinary);
# gid_io.WriteMesh((structure_model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);


for step in range(0, nsteps):

    time = Dt * step
    fluid_model_part.CloneTimeStep(time)
    structure_model_part.CloneTimeStep(time)

    print("time = ", time)

    # solving the fluid problem
    if(step < 1000):
        if(step > 3):
            print("solving only the fluid - starting procedure")
            fluid_solver.Solve()
            print("step ", step, " solved")
    else:
# if(step > 3):
        coupled_solver.Solve()
        PrintDisp(time, scalarout, tip)

    # print the results
# gid_io.WriteNodalResults(DISPLACEMENT,structure_model_part.Nodes,time,0)
# gid_io.WriteNodalResults(POSITIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)

    if(out == output_step):
# gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISPLACEMENT, fluid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
#        gid_io.WriteNodalResults(IS_INTERFACE,fluid_model_part.Nodes,time,0)

# gid_io.WriteNodalResults(IS_INTERFACE,structure_model_part.Nodes,time,0)
# gid_io.WriteNodalResults(POSITIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)
# gid_io.WriteNodalResults(NEGATIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)
# gid_io.WriteNodalResults(DISPLACEMENT,structure_model_part.Nodes,time,0)

        gid_io.Flush();
        out = 0

    out = out + 1
