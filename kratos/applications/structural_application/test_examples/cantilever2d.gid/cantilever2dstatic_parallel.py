from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import the configuration data as read from the GiD
import Kratos_Structural_Application_var

# find neighbours if required


def FindNeighbours():
    if(Kratos_Structural_Application_var.FindNodalNeighbours == "True"):
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        nodal_neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)
        nodal_neighbour_search.Execute()
    if(Kratos_Structural_Application_var.FindElementalNeighbours == "True"):
        neighbour_calculator = FindElementalNeighboursProcess(model_part, 2, 10)
        neighbour_calculator.Execute()

# importing the rotational dofs degrees of freedom if necessary


def RotationalDofs():
    if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
        for node in model_part.Nodes:
            node.AddDof(ROTATION_X)
            node.AddDof(ROTATION_Y)
            node.AddDof(ROTATION_Z)


#
#

from time import *
print(ctime())
t0 = clock()

# including kratos path
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MKLSolversApplication import *


# setting the domain size for the problem to be solved
domain_size = Kratos_Structural_Application_var.domain_size


# defining a model part
model_part = ModelPart("StructurePart")
model_part.AddNodalSolutionStepVariable(FORCE)
if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
    model_part.AddNodalSolutionStepVariable(ROTATION)


import structural_solver_static_parallel as SolverType
SolverType.AddVariables(model_part)


# reading a model
name = Kratos_Structural_Application_var.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh((model_part).GetMesh())
gid_io.FinalizeMesh()


# find neighbours if required
FindNeighbours()

model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D())
print("Linear elastic model selected")

print(model_part)
print(model_part.Properties)

# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

# importing the rotational dofs degrees of freedom if necessary
RotationalDofs()


# importing the solver files
SolverType.AddDofs(model_part)
solver = SolverType.StaticStructuralSolver(model_part, domain_size)
solver.structure_linear_solver = ParallelMKLPardisoSolver()
# solver.structure_linear_solver =  ParallelSkylineLUFactorizationSolver()


CT = Kratos_Structural_Application_var.Convergence_Tolerance;
AT = Kratos_Structural_Application_var.Absolute_Tolerance;

if(Kratos_Structural_Application_var.Convergence_Criteria == "Displacement_Criteria"):
    solver.conv_criteria = DisplacementCriteria(CT, AT)
elif(Kratos_Structural_Application_var.Convergence_Criteria == "Residual_Criteria"):
    solver.conv_criteria = ResidualCriteria(CT, AT)
elif(Kratos_Structural_Application_var.Convergence_Criteria == "And_Criteria"):
    Displacement = DisplacementCriteria(CT, AT)
    Residual = ResidualCriteria(CT, AT)
    solver.conv_criteria = AndCriteria(Residual, Displacement)
elif(Kratos_Structural_Application_var.Convergence_Criteria == "Or_Criteria"):
    Displacement = DisplacementCriteria(CT, AT)
    Residual = ResidualCriteria(CT, AT)
    solver.conv_criteria = OrCriteria(Residual, Displacement)


solver.Initialize()
solver.SetEchoLevel(1)


Dt = 0.02
nsteps = 2

mesh_name = 0.0
gid_io.InitializeMesh(mesh_name);
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()


print("initializing results")
gid_io.InitializeResults(mesh_name, (model_part).GetMesh())

for step in range(0, nsteps):
    time = Dt * step
    model_part.CloneTimeStep(time)
    if(step > 0):
        solver.Solve()
        # print the results
        gid_io.WriteNodalResults(DISPLACEMENT, model_part.Nodes, time, 0)
        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR, model_part, time)

gid_io.FinalizeResults()
tf = clock()
print("Analysis Completed. Time = ", tf - t0);
