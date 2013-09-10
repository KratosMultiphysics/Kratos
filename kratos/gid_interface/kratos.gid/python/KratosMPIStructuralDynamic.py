#
#
# import the configuration data as read from the GiD
import ProjectParameters


def PrintResults(model_part):
    print "Writing results. Please run Gid for viewing results of analysis."
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(
            varibles_dictionary[variable_name],
            model_part.Nodes,
            time,
            0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(
            varibles_dictionary[variable_name],
            model_part,
            time)


#
#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

#
#
# ATTENTION: here the order is important

from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *

if(ProjectParameters.LinearSolver == "SuperLUSolver"):
    from KratosMultiphysics.ExternalSolversApplication import *

if(ProjectParameters.SolverType == "ParallelSolver"):
    from KratosMultiphysics.MKLSolversApplication import *
# from now on the order is not anymore crucial
#
#

# defining variables to be used

varibles_dictionary = {"DISPLACEMENT": DISPLACEMENT,
                       "FORCE": FORCE,
                       "REACTION": REACTION,
                       "GREEN_LAGRANGE_STRAIN_TENSOR":
                       GREEN_LAGRANGE_STRAIN_TENSOR,
                       "ROTATION": ROTATION,
                       "PK2_STRESS_TENSOR": PK2_STRESS_TENSOR,
                       "MOMENT": MOMENT}

# defining a model part
model_part = ModelPart("StructurePart")
model_part.AddNodalSolutionStepVariable(FORCE)
if(ProjectParameters.Rotational_Dofs == "True"):
    model_part.AddNodalSolutionStepVariable(ROTATION)


# adding of Variables to Model Part should be here when the "very fix
# container will be ready"
if(ProjectParameters.Solution_method == "Newton-Raphson"):
    if(ProjectParameters.SolverType == "DynamicSolver"):
        if(ProjectParameters.Rotational_Dofs == "False"):
            import structural_solver_dynamic as SolverType
        else:
            import structural_solver_dynamic_rotation as SolverType


if(ProjectParameters.Solution_method == "LineSearch"):
    if(ProjectParameters.SolverType == "DynamicSolver"):
        if(ProjectParameters.Rotational_Dofs == "False"):
            import structural_solver_dynamic_general as SolverType
        else:
            import structural_solver_dynamic_rotation_general as SolverType


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
    nodal_neighbour_search = FindNodalNeighboursProcess(
        model_part, number_of_avg_elems, number_of_avg_nodes)
    nodal_neighbour_search.Execute()
if(ProjectParameters.FindElementalNeighbours == "True"):
    neighbour_calculator = FindElementalNeighboursProcess(model_part, 2, 10)
    neighbour_calculator.Execute()


print model_part
print model_part.Properties

# the buffer size should be set up here after the mesh is read for the
# first time
model_part.SetBufferSize(2)

# importing the rotational dofs degrees of freedom if necessary
if(ProjectParameters.Rotational_Dofs == "True"):
    for node in model_part.Nodes:
        node.AddDof(ROTATION_X)
        node.AddDof(ROTATION_Y)
        node.AddDof(ROTATION_Z)


SolverType.AddDofs(model_part)
solver = SolverType.DynamicStructuralSolver(model_part, domain_size)


# choosing the default value for the constitutive law
from materials import *
AssignMaterial(model_part.Properties)
print model_part.Properties

print "Linear elastic model selected"

#solver.structure_linear_solver = ProjectParameters.problem_name.LinearSolver()
if(ProjectParameters.LinearSolver == "SkylineLUFactorization"):
    solver.structure_linear_solver = SkylineLUFactorizationSolver()
elif(ProjectParameters.LinearSolver == "SuperLUSolver"):
    solver.structure_linear_solver = SuperLUSolver()
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
elif(ProjectParameters.LinearSolver == "ParallelMKLPardisoSolver"):
    solver.structure_linear_solver = ParallelMKLPardisoSolver()

CT = ProjectParameters.Convergence_Tolerance
AT = ProjectParameters.Absolute_Tolerance

if(ProjectParameters.Convergence_Criteria == "Displacement_Criteria"):
    solver.conv_criteria = DisplacementCriteria(CT, AT)
elif(ProjectParameters.Convergence_Criteria == "Residual_Criteria"):
    solver.conv_criteria = ResidualCriteria(CT, AT)
elif(ProjectParameters.Convergence_Criteria == "And_Criteria"):
    Displacement = DisplacementCriteria(CT, AT)
    Residual = ResidualCriteria(CT, AT)
    solver.conv_criteria = AndCriteria(Residual, Displacement)
elif(ProjectParameters.Convergence_Criteria == "Or_Criteria"):
    Displacement = DisplacementCriteria(CT, AT)
    Residual = ResidualCriteria(CT, AT)
    solver.conv_criteria = OrCriteria(Residual, Displacement)


solver.Initialize()
(solver.solver).SetEchoLevel(2)


Dt = ProjectParameters.Dt
MaxTime = ProjectParameters.max_time
Nsteps = ProjectParameters.nsteps
solver.max_iter = ProjectParameters.Max_Iter


gid_io.InitializeResults(mesh_name, (model_part).GetMesh())

for step in range(1, Nsteps):

    time = Dt * step
    model_part.CloneTimeStep(time)
    model_part.ProcessInfo[TIME_STEPS] = step

    print "STEP = ", step
    print "TIME = ", time
    # print model_part.ProcessInfo()[TIME]

    # solving the fluid problem
    if(step > 0):
        solver.Solve()
        PrintResults(model_part)

print "Analysis Completed "
gid_io.FinalizeResults()
