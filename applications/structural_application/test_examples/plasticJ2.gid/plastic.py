from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import the configuration data as read from the GiD
import Kratos_Structural_Application_var

#
#
# setting the domain size for the problem to be solved
domain_size = Kratos_Structural_Application_var.domain_size

#
#
# ATTENTION: here the order is important

# including kratos path
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *

if(Kratos_Structural_Application_var.LinearSolver == "SuperLUSolver"):
    from KratosMultiphysics.ExternalSolversApplication import *

# if(Kratos_Structural_Application_var.SolverType == "ParallelSolver"):
#     from KratosMultiphysics.MKLSolversApplication import *


# from now on the order is not anymore crucial
#
#

# defining a model part
model_part = ModelPart("StructurePart")
model_part.AddNodalSolutionStepVariable(FORCE)
if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
    model_part.AddNodalSolutionStepVariable(ROTATION)


# adding of Variables to Model Part should be here when the "very fix container will be ready"
if(Kratos_Structural_Application_var.SolverType == "StaticSolver"):
    import structural_solver_static
    structural_solver_static.AddVariables(model_part)
elif(Kratos_Structural_Application_var.SolverType == "DynamicSolver"):
    if(Kratos_Structural_Application_var.Rotational_Dofs == "False"):
        import structural_solver_dynamic
        structural_solver_dynamic.AddVariables(model_part)
    else:
        import structural_solver_dynamic_rotation
        structural_solver_dynamic_rotation.AddVariables(model_part)
elif(Kratos_Structural_Application_var.SolverType == "ParallelSolver"):
    import structural_solver_static_parallel
    structural_solver_static_parallel.AddVariables(model_part)
elif(Kratos_Structural_Application_var.SolverType == "ArcLengthSolver"):
    import structural_solver_static_arc_length
    structural_solver_static_arc_length.AddVariables(model_part)
elif(Kratos_Structural_Application_var.SolverType == "LineSearchesSolver"):
    import structural_solver_static_general
    structural_solver_static_general.AddVariables(model_part)

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
if(Kratos_Structural_Application_var.FindNodalNeighbours == "True"):
    number_of_avg_elems = 10
    number_of_avg_nodes = 10
    nodal_neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)
    nodal_neighbour_search.Execute()
if(Kratos_Structural_Application_var.FindElementalNeighbours == "True"):
    neighbour_calculator = FindElementalNeighboursProcess(model_part, 2, 10)
    neighbour_calculator.Execute()


print(model_part)
print(model_part.Properties)

# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

# importing the rotational dofs degrees of freedom if necessary
if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
    for node in model_part.Nodes:
        node.AddDof(ROTATION_X)
        node.AddDof(ROTATION_Y)
        node.AddDof(ROTATION_Z)

# importing the solver files
if(Kratos_Structural_Application_var.SolverType == "StaticSolver"):
    structural_solver_static.AddDofs(model_part)
    solver = structural_solver_static.StaticStructuralSolver(model_part, domain_size)
elif(Kratos_Structural_Application_var.SolverType == "DynamicSolver"):
    if(Kratos_Structural_Application_var.Rotational_Dofs == "False"):
        structural_solver_dynamic.AddDofs(model_part)
        solver = structural_solver_dynamic.DynamicStructuralSolver(model_part, domain_size)
    else:
        structural_solver_dynamic_rotation.AddDofs(model_part)
        solver = structural_solver_dynamic_rotation.DynamicStructuralSolver(model_part, domain_size)
elif(Kratos_Structural_Application_var.SolverType == "ParallelSolver"):
    structural_solver_static_parallel.AddDofs(model_part)
    solver = structural_solver_static_parallel.StaticStructuralSolver(model_part, domain_size)
elif(Kratos_Structural_Application_var.SolverType == "ArcLengthSolver"):
    structural_solver_static_arc_length.AddDofs(model_part)
    solver = structural_solver_static_arc_length.StaticStructuralSolver(model_part, domain_size)
    model_part.ProcessInfo[LAMNDA] = 0.00
elif(Kratos_Structural_Application_var.SolverType == "LineSearchesSolver"):
    structural_solver_static_general.AddDofs(model_part)
    solver = structural_solver_static_general.StaticStructuralSolver(model_part, domain_size)

# choosing the default value for the constitutive law
if(domain_size == 2):
    for prop in model_part.Properties:
        prop.SetValue(CONSTITUTIVE_LAW, PlaneStressJ2())  # Isotropic2D() )#PlaneStressJ2()
        prop.SetValue(DENSITY, 2700.0000);
        prop.SetValue(POISSON_RATIO, 0.3);
        prop.SetValue(YIELD_STRESS, 2.76e8);
        prop.SetValue(PLASTIC_MODULUS, 0.0);  # 6.34317e8);
        prop.SetValue(YOUNG_MODULUS, 70000e6);

else:
    for prop in model_part.Properties:
        prop.SetValue(CONSTITUTIVE_LAW, Isotropic3D())
# creating a fluid solver object
# model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
# model_part.Properties[2].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )

print("Linear elastic model selected")

# solver.structure_linear_solver = Kratos_Structural_Application_var.problem_name.LinearSolver()
if(Kratos_Structural_Application_var.LinearSolver == "SkylineLUFactorization"):
    solver.structure_linear_solver = SkylineLUFactorizationSolver()
elif(Kratos_Structural_Application_var.LinearSolver == "SuperLUSolver"):
    solver.structure_linear_solver = SuperLUSolver()
elif(Kratos_Structural_Application_var.LinearSolver == "BiCGStab_ILU0"):
    pILUPrecond = ILU0Preconditioner()
    LST = Kratos_Structural_Application_var.Linear_Solver_Tolerance
    LSMI = Kratos_Structural_Application_var.Linear_Solver_Max_Iteration
    solver.structure_linear_solver = BICGSTABSolver(LST, LSMI, pILUPrecond)
elif(Kratos_Structural_Application_var.LinearSolver == "BiCGStab_DIAG"):
    pDiagPrecond = DiagonalPreconditioner()
    LST = Kratos_Structural_Application_var.Linear_Solver_Tolerance
    LSMI = Kratos_Structural_Application_var.Linear_Solver_Max_Iteration
    solver.structure_linear_solver = BICGSTABSolver(LST, LSMI, pDiagPrecond)
# elif(Kratos_Structural_Application_var.LinearSolver == "ParallelMKLPardisoSolver"):
#     solver.structure_linear_solver = ParallelMKLPardisoSolver()

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
(solver.solver).SetEchoLevel(2);


Dt = Kratos_Structural_Application_var.Dt
MaxTime = Kratos_Structural_Application_var.max_time
Nsteps = Kratos_Structural_Application_var.nsteps
solver.max_iter = Kratos_Structural_Application_var.Max_Iter

for node in model_part.Nodes:
    node.SetSolutionStepValue(FORCE_X, 0, 0.0)
    node.SetSolutionStepValue(FORCE_Y, 0, 0.0000)  # 7217.5
    node.SetSolutionStepValue(FORCE_Z, 0, 0.0)

for node in model_part.Nodes:
    if(node.X == 7.0):
        node.SetSolutionStepValue(DISPLACEMENT_X, 0, 0.0132)
        node.Fix(DISPLACEMENT_X)
        node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Y)


gid_io.InitializeResults(mesh_name, (model_part).GetMesh())

for step in range(1, Nsteps):

    time = Dt * step
    model_part.CloneTimeStep(time)
    model_part.ProcessInfo[TIME_STEPS] = step

    print("STEP = ", step)
    print("TIME = ", time)
    # print model_part.ProcessInfo()[TIME]

    # solving the fluid problem
    if(step > 0):
        solver.Solve()
        if(Kratos_Structural_Application_var.SolverType == "ArcLengthSolver"):
            structural_solver_static_arc_length.ChangeCondition(model_part, model_part.ProcessInfo[LAMNDA])
            print(model_part.ProcessInfo[LAMNDA]);

        print("Writing results. Please run Gid for viewing results of analysis.")
        # gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
        # gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time,domain_size)
        gid_io.WriteNodalResults(FORCE, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(REACTION, model_part.Nodes, time, 0)
        # gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR,model_part,time)
        gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, model_part, time)
        # gid_io.PrintOnGaussPoints(DAMAGE,model_part,time)
        gid_io.WriteNodalResults(DISPLACEMENT, model_part.Nodes, time, 0)
        if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
            gid_io.WriteNodalResults(ROTATION, model_part.Nodes, time, 0)
        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR, model_part, time)
        if(Kratos_Structural_Application_var.SolverType == "DynamicSolver"):
            gid_io.WriteNodalResults(VELOCITY, model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(ACCELERATION, model_part.Nodes, time, 0)

print("Analysis Completed ")
gid_io.FinalizeResults()
