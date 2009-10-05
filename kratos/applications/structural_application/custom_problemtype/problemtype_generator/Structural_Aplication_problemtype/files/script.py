#import the configuration data as read from the GiD
import Structural_Aplication_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = Structural_Aplication_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = Structural_Aplication_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = Structural_Aplication_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_StructuralApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *

if(Structural_Aplication_var.LinearSolver == "SuperLUSolver"):
    from KratosExternalSolversApplication import *

## from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part
model_part = ModelPart("StructurePart");  
model_part.AddNodalSolutionStepVariable(FORCE);
if(Structural_Aplication_var.Rotational_Dofs == "True"):
  model_part.AddNodalSolutionStepVariable(ROTATION);


#adding of Variables to Model Part should be here when the "very fix container will be ready"
if(Structural_Aplication_var.SolverType == "StaticSolver"): 
    import structural_solver_static
    structural_solver_static.AddVariables(model_part)
elif(Structural_Aplication_var.SolverType == "DynamicSolver"):
    if(Structural_Aplication_var.Rotational_Dofs == "False"):
	import structural_solver_dynamic
	structural_solver_dynamic.AddVariables(model_part)
    else:
	import structural_solver_dynamic_rotation
	structural_solver_dynamic_rotation.AddVariables(model_part)
elif(Structural_Aplication_var.SolverType == "ParallelSolver"):
    import structural_solver_static_parallel
    structural_solver_static_parallel.AddVariables(model_part)
elif(Structural_Aplication_var.SolverType == "ArcLengthSolver"):
    import structural_solver_static_arc_length
    structural_solver_static_arc_length.AddVariables(model_part)
elif(Structural_Aplication_var.SolverType == "LineSearchesSolver"):
    import structural_solver_static_general
    structural_solver_static_general.AddVariables(model_part)

#reading a model
name = Structural_Aplication_var.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

##find neighbours if required
if(Structural_Aplication_var.FindNodalNeighbours == "True"):
    number_of_avg_elems = 10
    number_of_avg_nodes = 10
    nodal_neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
    nodal_neighbour_search.Execute()
if(Structural_Aplication_var.FindElementalNeighbours == "True"):
    neighbour_calculator = FindElementalNeighboursProcess(model_part,2,10);
    neighbour_calculator.Execute()


print model_part
print model_part.Properties

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

##importing the rotational dofs degrees of freedom if necessary
if(Structural_Aplication_var.Rotational_Dofs == "True"):
  for node in model_part.Nodes:
    node.AddDof(ROTATION_X)
    node.AddDof(ROTATION_Y)
    node.AddDof(ROTATION_Z)

#importing the solver files
if(Structural_Aplication_var.SolverType == "StaticSolver"):
    structural_solver_static.AddDofs(model_part)
    solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size) 
elif(Structural_Aplication_var.SolverType == "DynamicSolver"):
    if(Structural_Aplication_var.Rotational_Dofs == "False"):
	structural_solver_dynamic.AddDofs(model_part)
	solver = structural_solver_dynamic.DynamicStructuralSolver(model_part,domain_size)
    else:
	structural_solver_dynamic_rotation.AddDofs(model_part)
	solver = structural_solver_dynamic_rotation.DynamicStructuralSolver(model_part,domain_size)
elif(Structural_Aplication_var.SolverType == "ParallelSolver"):
    structural_solver_static_parallel.AddDofs(model_part)
    solver = structural_solver_static_parallel.StaticStructuralSolver(model_part,domain_size)
elif(Structural_Aplication_var.SolverType == "ArcLengthSolver"):
    structural_solver_static_arc_length.AddDofs(model_part)
    solver = structural_solver_static_arc_length.StaticStructuralSolver(model_part,domain_size)
elif(Structural_Aplication_var.SolverType == "LineSearchesSolver"):
    structural_solver_static_general.AddDofs(model_part)
    solver = structural_solver_static_general.StaticStructuralSolver(model_part,domain_size)

##choosing the default value for the constitutive law 
if(domain_size == 2):
  for prop in model_part.Properties:
    prop.SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
else:
  for prop in model_part.Properties:
    prop.SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
#creating a fluid solver object
#model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
#model_part.Properties[2].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )

print "Linear elastic model selected"

#solver.structure_linear_solver = Structural_Aplication_var.problem_name.LinearSolver()  
if(Structural_Aplication_var.LinearSolver == "SkylineLUFactorization"):
    solver.structure_linear_solver  =  SkylineLUFactorizationSolver()
elif(Structural_Aplication_var.LinearSolver == "SuperLUSolver"):
    solver.structure_linear_solver  =  SuperLUSolver()
elif(Structural_Aplication_var.LinearSolver == "BiCGStab_ILU0"):
    pILUPrecond = ILU0Preconditioner() 
    LST  = Structural_Aplication_var.Linear_Solver_Tolerance
    LSMI = Structural_Aplication_var.Linear_Solver_Max_Iteration 
    solver.structure_linear_solver  =  BICGSTABSolver(LST, LSMI, pILUPrecond)
elif(Structural_Aplication_var.LinearSolver == "BiCGStab_DIAG"):
    pDiagPrecond = DiagonalPreconditioner()
    LST  = Structural_Aplication_var.Linear_Solver_Tolerance
    LSMI = Structural_Aplication_var.Linear_Solver_Max_Iteration  
    solver.structure_linear_solver  =  BICGSTABSolver(LST,LSMI,pDiagPrecond)



CT = Structural_Aplication_var.Convergence_Tolerance;
AT = Structural_Aplication_var.Absolute_Tolerance; 

if(Structural_Aplication_var.Convergence_Criteria == "Displacement_Criteria"):
    solver.conv_criteria  =  DisplacementCriteria(CT,AT)
elif(Structural_Aplication_var.Convergence_Criteria == "Residual_Criteria"): 
    solver.conv_criteria  =   ResidualCriteria(CT,AT)
elif(Structural_Aplication_var.Convergence_Criteria == " Both_Criteria"): 
    Displacement   =  DisplacementCriteria(CT,AT)
    Residual       =   ResidualCriteria(CT,AT)
    solver.conv_criteria  = ResDisCriteria(Residual, Displacement)
        

  
solver.Initialize()
(solver.solver).SetEchoLevel(2);


Dt      = Structural_Aplication_var.Dt
MaxTime = Structural_Aplication_var.max_time
Nsteps  = Structural_Aplication_var.nsteps
solver.max_iter =  Structural_Aplication_var.Max_Iter


gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

for step in range(0,Nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    print "STEP = ", step
    print "TIME = ", time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        solver.Solve()

        print "Writing results. Please run Gid for viewing results of analysis."
        #gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
        #gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time,domain_size)
	gid_io.WriteNodalResults(FORCE,model_part.Nodes,time,0)
	gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
	gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR,model_part,time)
	#gid_io.PrintOnGaussPoints(DAMAGE,model_part,time)
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
	if(Structural_Aplication_var.Rotational_Dofs == "True"):
	  gid_io.WriteNodalResults(ROTATION,model_part.Nodes,time,0)
        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
        if(Structural_Aplication_var.SolverType == "DynamicSolver"):
	  gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
	  gid_io.WriteNodalResults(ACCELERATION,model_part.Nodes,time,0)

print "Analysis Completed "
gid_io.FinalizeResults()

          
        

