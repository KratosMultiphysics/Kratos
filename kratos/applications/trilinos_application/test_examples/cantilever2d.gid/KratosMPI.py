#importing MPI ... for this boost 1.35 or superior is needed
import mpi

##################################################################
##################################################################
#import the configuration data as read from the GiD
import ProjectParameters



def PrintResults(model_part):
        print "Writing results. Please run Gid for viewing results of analysis."
        for variable_name in ProjectParameters.nodal_results:
            gid_io.WriteNodalResults(varibles_dictionary[variable_name],model_part.Nodes,time,0)
        for variable_name in ProjectParameters.gauss_points_results:
            gid_io.PrintOnGaussPoints(varibles_dictionary[variable_name],model_part,time)



##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size


##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = ProjectParameters.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = ProjectParameters.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_StructuralApplication = True
applications_interface.Import_KratosMetisApplication = True
applications_interface.Import_KratosTrilinosApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *
from KratosMetisApplication import *
from KratosTrilinosApplication import *

if(ProjectParameters.LinearSolver == "SuperLUSolver"):
    from KratosExternalSolversApplication import *

if(ProjectParameters.SolverType == "ParallelSolver"):
    applications_interface.Import_KratosMKLSolversApplication = True
    from KratosMKLSolversApplication import *
## from now on the order is not anymore crucial
##################################################################
##################################################################

## defining variables to be used

varibles_dictionary = {"DISPLACEMENT" : DISPLACEMENT,
                       "FORCE" : FORCE,
                       "REACTION" : REACTION,
                       "GREEN_LAGRANGE_STRAIN_TENSOR" : GREEN_LAGRANGE_STRAIN_TENSOR,
                       "ROTATION" : ROTATION,
                       "PK2_STRESS_TENSOR" : PK2_STRESS_TENSOR}

#defining a model part
model_part = ModelPart("StructurePart");  


model_part.AddNodalSolutionStepVariable(FORCE);
if(ProjectParameters.Rotational_Dofs == "True"):
  model_part.AddNodalSolutionStepVariable(ROTATION);


#adding of Variables to Model Part should be here when the "very fix container will be ready"
if(ProjectParameters.SolverType == "StaticSolver"): 
    import trilinos_structural_solver_static
    trilinos_structural_solver_static.AddVariables(model_part)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
elif(ProjectParameters.SolverType == "ParallelSolver"):
    import structural_solver_static_parallel
    structural_solver_static_parallel.AddVariables(model_part)
elif(ProjectParameters.SolverType == "ArcLengthSolver"):
    import structural_solver_static_arc_length
    structural_solver_static_arc_length.AddVariables(model_part)
elif(ProjectParameters.SolverType == "LineSearchesSolver"):
    import structural_solver_static_general
    structural_solver_static_general.AddVariables(model_part)

#reading a model
name = ProjectParameters.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)

number_of_partitions = mpi.size #we set it equal to the number of processors
if( mpi.rank == 0):
	print "number_of_partitions", number_of_partitions
	
partitioner = MetisPartitioningProcess(model_part, model_part_io, number_of_partitions, domain_size);
partitioner.Execute()

mesh_name = mpi.rank
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

##find neighbours if required
if(ProjectParameters.FindNodalNeighbours == "True"):
    number_of_avg_elems = 10
    number_of_avg_nodes = 10
    nodal_neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
    nodal_neighbour_search.Execute()
if(ProjectParameters.FindElementalNeighbours == "True"):
    neighbour_calculator = FindElementalNeighboursProcess(model_part,2,10);
    neighbour_calculator.Execute()


print model_part
print model_part.Properties

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

##importing the rotational dofs degrees of freedom if necessary
if(ProjectParameters.Rotational_Dofs == "True"):
  for node in model_part.Nodes:
    node.AddDof(ROTATION_X)
    node.AddDof(ROTATION_Y)
    node.AddDof(ROTATION_Z)

#importing the solver files
if(ProjectParameters.SolverType == "StaticSolver"):
    trilinos_structural_solver_static.AddDofs(model_part)
    solver = trilinos_structural_solver_static.StaticStructuralSolver(model_part,domain_size) 
elif(ProjectParameters.SolverType == "ParallelSolver"):
    structural_solver_static_parallel.AddDofs(model_part)
    solver = structural_solver_static_parallel.StaticStructuralSolver(model_part,domain_size)
elif(ProjectParameters.SolverType == "ArcLengthSolver"):
    structural_solver_static_arc_length.AddDofs(model_part)
    solver = structural_solver_static_arc_length.StaticStructuralSolver(model_part,domain_size)
    model_part.ProcessInfo[LAMNDA] = 0.00;  
elif(ProjectParameters.SolverType == "LineSearchesSolver"):
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
########################################################
#defining the linear solver
aztec_parameters = ParameterList()
aztec_parameters.set("AZ_solver","AZ_gmres");
aztec_parameters.set("AZ_kspace",100);
aztec_parameters.set("AZ_output",32);

preconditioner_type = "Amesos"
preconditioner_parameters = ParameterList()
preconditioner_parameters.set("amesos: solver type", "Amesos_Klu");

##preconditioner_type = "ILU"
##preconditioner_parameters = ParameterList()

overlap_level = 3
nit_max = 300
tol = 1e-6

solver.structure_linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,tol,nit_max,overlap_level);


CT = ProjectParameters.Convergence_Tolerance;
AT = ProjectParameters.Absolute_Tolerance; 

if(ProjectParameters.Convergence_Criteria == "Displacement_Criteria"):
    solver.conv_criteria  =  TrilinosDisplacementCriteria(CT,AT, solver.Comm)
elif(ProjectParameters.Convergence_Criteria == "Residual_Criteria"): 
    solver.conv_criteria  =   TrilinosResidualCriteria(CT,AT)
elif(ProjectParameters.Convergence_Criteria == "And_Criteria"): 
    Displacement   =   TrilinosDisplacementCriteria(CT,AT)
    Residual       =   TrilinosResidualCriteria(CT,AT)
    solver.conv_criteria  = TrilinosAndCriteria(Residual, Displacement)
elif(ProjectParameters.Convergence_Criteria == "Or_Criteria"): 
    Displacement   =   TrilinosDisplacementCriteria(CT,AT)
    Residual       =   TrilinosResidualCriteria(CT,AT)
    solver.conv_criteria  = TrilinosOrCriteria(Residual, Displacement)
        


  
solver.Initialize()
(solver.solver).SetEchoLevel(2);

time = 1.00;
step = 1;

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

model_part.CloneTimeStep(time)
model_part.ProcessInfo[TIME_STEPS] = step

solver.Solve()
if(ProjectParameters.SolverType == "ArcLengthSolver"):
	structural_solver_static_arc_length.ChangeCondition(model_part, model_part.ProcessInfo[LAMNDA])
	print model_part.ProcessInfo[LAMNDA];

PrintResults(model_part)

print "Analysis Completed for process #", mpi.rank
gid_io.FinalizeResults()

          
        

