from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters


def ChangeCondition(model_part, lamda):
    for node in model_part.Nodes:
        new_load = node.GetSolutionStepValue(POINT_LOAD) * lamda
        node.SetSolutionStepValue(POINT_LOAD, 0, new_load)


def PrintResults(model_part):
    print("Writing results. Please run Gid for viewing results of analysis.")
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(varibles_dictionary[variable_name], model_part.Nodes, time, 0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(varibles_dictionary[variable_name], model_part, time)


#
#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

#
#
# ATTENTION: here the order is important


from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *

if(ProjectParameters.LinearSolver == "SuperLUSolver" or ProjectParameters.LinearSolver == "SuperLUIterativeSolver"):
    from KratosMultiphysics.ExternalSolversApplication import *

if(ProjectParameters.LinearSolver == "ParallelMKLPardisoSolver" or ProjectParameters.LinearSolver == "MKLPardisoSolver"):
    from KratosMultiphysics.MKLSolversApplication import *


# from now on the order is not anymore crucial
#
#

# defining variables to be used

varibles_dictionary = {"DISPLACEMENT": DISPLACEMENT,
                       "POINT_LOAD": POINT_LOAD,
                       "REACTION": REACTION,
                       "GREEN_LAGRANGE_STRAIN_TENSOR": GREEN_LAGRANGE_STRAIN_TENSOR,
                       "ROTATION": ROTATION,
                       "PK2_STRESS_TENSOR": PK2_STRESS_TENSOR,
                       "MOMENT": MOMENT}

# defining a model part
model_part = ModelPart("StructurePart")
model_part.AddNodalSolutionStepVariable(POINT_LOAD)
if(ProjectParameters.Rotational_Dofs == "True"):
    model_part.AddNodalSolutionStepVariable(ROTATION)
    model_part.AddNodalSolutionStepVariable(MOMENTUM)


# adding of Variables to Model Part should be here when the "very fix container will be ready"

if(ProjectParameters.Solution_method == "Newton-Raphson"):
    if(ProjectParameters.SolverType == "StaticSolver"):
        if(ProjectParameters.LinearSolver == "ParallelMKLPardisoSolver"):
            import structural_solver_static_parallel as SolverType
        else:
            import structural_solver_static as SolverType

if(ProjectParameters.Solution_method == "ArcLength"):
    import structural_solver_static_arc_length as SolverType

if(ProjectParameters.Solution_method == "LineSearch"):
    import structural_solver_static_general as SolverType

SolverType.AddVariables(model_part)

# reading a model
name = ProjectParameters.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
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
if(ProjectParameters.FindNodalNeighbours == "True"):
    number_of_avg_elems = 10
    number_of_avg_nodes = 10
    nodal_neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)
    nodal_neighbour_search.Execute()
if(ProjectParameters.FindElementalNeighbours == "True"):
    neighbour_calculator = FindElementalNeighboursProcess(model_part, 2, 10);
    neighbour_calculator.Execute()


print(model_part)

# choosing the default value for the constitutive law
from materials import *
AssignMaterial(model_part.Properties)
print(model_part.Properties)


# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

# importing the rotational dofs degrees of freedom if necessary
if(ProjectParameters.Rotational_Dofs == "True"):
    for node in model_part.Nodes:
        node.AddDof(ROTATION_X, MOMENTUM_X);
        node.AddDof(ROTATION_Y, MOMENTUM_Y);
        node.AddDof(ROTATION_Z, MOMENTUM_Z);


# importing the solver files
SolverType.AddDofs(model_part)
solver = SolverType.StaticStructuralSolver(model_part, domain_size)

if(ProjectParameters.Solution_method == "ArcLength"):
    model_part.ProcessInfo[LAMNDA] = 0.00;


# solver.structure_linear_solver = ProjectParameters.problem_name.LinearSolver()
if(ProjectParameters.LinearSolver == "SkylineLUFactorization"):
    solver.structure_linear_solver = SkylineLUFactorizationSolver()
elif(ProjectParameters.LinearSolver == "SuperLUSolver"):
    solver.structure_linear_solver = SuperLUSolver()
elif(ProjectParameters.LinearSolver == "SuperLUIterativeSolver"):
    solver.structure_linear_solver = SuperLUIterativeSolver()
elif(ProjectParameters.LinearSolver == "BiCGStab_ILU0"):
    pILUPrecond = ILU0Preconditioner()
    LST = ProjectParameters.Linear_Solver_Tolerance
    LSMI = ProjectParameters.Linear_Solver_Max_Iteration
    solver.structure_linear_solver = BICGSTABSolver(LST, LSMI, pILUPrecond)
elif(ProjectParameters.LinearSolver == "BiCGStab_DIAG"):
    pDiagPrecond = DiagonalPreconditioner()
    LST = ProjectParameters.Linear_Solver_Tolerance
    LSMI = ProjectParameters.Linear_Solver_Max_Iteration
    solver.structure_linear_solver = BICGSTABSolver(LST, LSMI, pDiagPrecond)
elif(ProjectParameters.LinearSolver == "MKLPardisoSolver"):
    solver.structure_linear_solver = MKLPardisoSolver()
elif(ProjectParameters.LinearSolver == "ParallelMKLPardisoSolver"):
    solver.structure_linear_solver = ParallelMKLPardisoSolver()


print(solver.structure_linear_solver)

DCT = ProjectParameters.Displacement_Convergence_Tolerance;
DAT = ProjectParameters.Displacement_Absolute_Tolerance;
RCT = ProjectParameters.Residual_Convergence_Tolerance;
RAT = ProjectParameters.Residual_Absolute_Tolerance;


if(ProjectParameters.Convergence_Criteria == "Displacement_Criteria"):
    solver.conv_criteria = DisplacementCriteria(DCT, DAT)
elif(ProjectParameters.Convergence_Criteria == "Residual_Criteria"):
    solver.conv_criteria = ResidualCriteria(RCT, RAT)
elif(ProjectParameters.Convergence_Criteria == "And_Criteria"):
    Displacement = DisplacementCriteria(DCT, DAT)
    Residual = ResidualCriteria(RCT, RAT)
    solver.conv_criteria = AndCriteria(Residual, Displacement)
elif(ProjectParameters.Convergence_Criteria == "Or_Criteria"):
    Displacement = DisplacementCriteria(DCT, DAT)
    Residual = ResidualCriteria(RCT, RAT)
    solver.conv_criteria = OrCriteria(Residual, Displacement)


for node in model_part.Nodes:
    if(node.IsFixed(DISPLACEMENT_X)):
        if(node.X > -0.0001 and node.Z < 0.0001 and node.Z < 25.0):
            node.Fix(DISPLACEMENT_Z)
        else:
            node.Free(DISPLACEMENT_Z)


# setting to ensure a linear elastic solution
solver.MaxNewtonRapshonIterations = 1
solver.MoveMeshFlag = False
solver.Initialize()

(solver.solver).SetEchoLevel(2);

time = 1.00;
step = 1;


def FindMinDispY(nodes):
    dx_min = 0.0
    dy_min = 0.0
    dz_min = 0.0
    for node in nodes:
        dx = node.GetSolutionStepValue(DISPLACEMENT_X)
        dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
        dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
        if(dx_min > dx):
            dx_min = dx
        if(dy_min > dy):
            dy_min = dy
        if(dz_min > dz):
            dz_min = dz

    return [dx_min, dy_min, dz_min]


gid_io.InitializeResults(mesh_name, (model_part).GetMesh())
print("Result files initialized ")
model_part.CloneTimeStep(time)
model_part.ProcessInfo[TIME_STEPS] = step

solver.Solve()

[dx_min, dy_min, dz_min] = FindMinDispY(model_part.Nodes)
print([dx_min, dy_min, dz_min])

import sys
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking

if (benchmarking.InBenchmarkingMode()):
    abs_tol = 1e-9
    rel_tol = 1e-5
    benchmarking.Output(dx_min, "min dx", abs_tol, rel_tol)
    benchmarking.Output(dy_min, "min dy", abs_tol, rel_tol)
    benchmarking.Output(dz_min, "min dz", abs_tol, rel_tol)


if(ProjectParameters.Solution_method == "ArcLength"):
    ChangeCondition(model_part, model_part.ProcessInfo[LAMNDA])
    print(model_part.ProcessInfo[LAMNDA]);

PrintResults(model_part)

print("Analysis Completed ")
gid_io.FinalizeResults()
