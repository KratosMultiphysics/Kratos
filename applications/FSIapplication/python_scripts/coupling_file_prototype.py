from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# setting the domain size for the problem to be solved
domain_size = 3


from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.ALEApplication import *


# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
structure_model_part = ModelPart("StructurePart")

# introducing input file name
input_file_name = "file_name "

# reading the fluid part
gid_io = GidIO(input_file_name + str("_fluid"), GiDPostMode.GiD_PostBinary)
gid_io.ReadMesh(fluid_model_part.GetMesh())
print(fluid_model_part)
print("fluid model read correctly")

# reading the structural part
data_io = DatafileIO(input_file_name + str("_structure"))
data_io.ReadMesh(structure_model_part.GetMesh())
print(structure_model_part)
print("structural model read correctly")

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)
structure_model_part.SetBufferSize(2)

# assigning the fluid properties
viscosity = 0.00001
density = 1.21
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY, 0, viscosity)
    node.SetSolutionStepValue(DENSITY, 0, density)

#
# importing the solvers needed
import mesh_solver
import incompressible_fluid_solver
import NonConformant_OneSideMap
# import Conformant_OneSideMap
import structural_solver_dynamic
import ExplicitCoupling

incompressible_fluid_solver.AddVariables(fluid_model_part)
mesh_solver.AddVariables(fluid_model_part)
NonConformant_OneSideMap.AddVariables(fluid_model_part, structure_model_part)
structural_solver_dynamic.AddVariables(structure_model_part)
print("solution variables added correctly")

# creating the solvers
# fluid solver
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part, domain_size)
fluid_solver.laplacian_form = 1
# standard laplacian form
fluid_solver.Initialize()
print("fluid solver created")

# mesh solver
reform_dofs_at_each_step = False
mesh_solver = mesh_solver.MeshSolver(fluid_model_part, domain_size, reform_dofs_at_each_step)
mesh_solver.Initialize()
mesh_solver.solver.SetEchoLevel(0)
print("mesh solver created")

# structure solver
structure_solver = structural_solver_dynamic.DynamicStructuralSolver(structure_model_part, domain_size)
structure_solver.echo_level = 0
pILUPrecond = ILU0Preconditioner()
structure_solver.structure_linear_solver = BICGSTABSolver(1e-8, 5000, pILUPrecond)
structure_solver.Initialize()
print("structural solver created")

# mapper
# non conformant mapper
mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(fluid_model_part, structure_model_part)
# conformant point to point
# utilities = VariableUtils()
# interface_fluid_nodes = (utilities).SelectNodeList(IS_INTERFACE,1.0,fluid_model_part.Nodes)
# interface_structure_nodes = (utilities).SelectNodeList(IS_INTERFACE,1.0,structure_model_part.Nodes)
# mapper = Conformant_OneSideMap.Conformant_OneSideMap(interface_fluid_nodes,interface_structure_nodes)

print("mapper created")

# creating the coupled solver
coupled_solver = ExplicitCoupling.ExplicitCoupling(fluid_model_part, structure_model_part, structure_solver, fluid_solver, mesh_solver, mapper, domain_size)
print("coupled solver created")


# settings to be changed
nsteps = 2000
output_step = 1

Dt = 0.01

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
    if(step < 100):
        if(step > 3):
            print("solving only the fluid - starting procedure")
            fluid_solver.Solve()
    else:
        coupled_solver.Solve()

    # print the results
# gid_io.WriteNodalResults(DISPLACEMENT,structure_model_part.Nodes,time,0)
# gid_io.WriteNodalResults(POSITIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)

# gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time,0)
    gid_io.WriteNodalResults(VELOCITY, fluid_model_part.Nodes, time, 0)
    gid_io.WriteNodalResults(DISPLACEMENT, fluid_model_part.Nodes, time, 0)
#    gid_io.WriteNodalResults(POSITIVE_FACE_PRESSURE,structure_model_part.Nodes,time,0)
    gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)

    gid_io.Flush();
