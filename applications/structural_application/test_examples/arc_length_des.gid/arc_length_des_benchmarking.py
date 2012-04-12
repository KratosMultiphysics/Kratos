
def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node
    
def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(REACTION_Y), "Node 1 Reaction_y", 1.0)
    

#import the configuration data as read from the GiD
import Kratos_Structural_Application_var

##find neighbours if required
def FindNeighbours():
  if(Kratos_Structural_Application_var.FindNodalNeighbours == "True"):
    number_of_avg_elems = 10
    number_of_avg_nodes = 10
    nodal_neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
    nodal_neighbour_search.Execute()
  if(Kratos_Structural_Application_var.FindElementalNeighbours == "True"):
    neighbour_calculator = FindElementalNeighboursProcess(model_part,2,10);
    neighbour_calculator.Execute()
  
##importing the rotational dofs degrees of freedom if necessary
def RotationalDofs(): 
  if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
    for node in model_part.Nodes:
      node.AddDof(ROTATION_X)
      node.AddDof(ROTATION_Y)
      node.AddDof(ROTATION_Z)


def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node


##################################################################
##################################################################
from time import *
print ctime()
t0 = clock()

#including kratos path
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
if(Kratos_Structural_Application_var.LinearSolver == "SuperLUSolver"):
    from KratosMultiphysics.ExternalSolversApplication import *
if(Kratos_Structural_Application_var.SolverType == "ParallelSolver"):
    from KratosMultiphysics.MKLSolversApplication import *
    
    
#setting the domain size for the problem to be solved
domain_size = Kratos_Structural_Application_var.domain_size


#defining a model part
model_part = ModelPart("StructurePart");  
model_part.AddNodalSolutionStepVariable(FORCE);
if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
  model_part.AddNodalSolutionStepVariable(ROTATION);


#adding of Variables to Model Part should be here when the "very fix container will be ready"
if(Kratos_Structural_Application_var.SolverType == "StaticSolver"): 
    import structural_solver_static as SolverType
elif(Kratos_Structural_Application_var.SolverType == "DynamicSolver"):
     if(Kratos_Structural_Application_var.TimeIntegration=="Bossak_Newmark"):
        if(Kratos_Structural_Application_var.Rotational_Dofs == "False"):
	      import structural_solver_dynamic as SolverType 
	else:
	     import structural_solver_dynamic_rotation as SolverType 
     elif(Kratos_Structural_Application_var.TimeIntegration=="Central_Differences"):
        import structural_solver_dynamic_central_differences as SolverType 
elif(Kratos_Structural_Application_var.SolverType == "ArcLengthSolver"):
    import structural_solver_static_arc_length as SolverType 
elif(Kratos_Structural_Application_var.SolverType == "LineSearchesSolver"):
    import structural_solver_static_general as SolverType
elif(Kratos_Structural_Application_var.SolverType == "ParallelSolver"):
    import structural_solver_static_parallel as SolverType     



SolverType.AddVariables(model_part)

#reading a model
name = Kratos_Structural_Application_var.problem_name

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
FindNeighbours();

print model_part
print model_part.Properties

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

##importing the rotational dofs degrees of freedom if necessary
RotationalDofs()



#importing the solver files
SolverType.AddDofs(model_part)

if(Kratos_Structural_Application_var.SolverType == "StaticSolver"):
    solver = SolverType.StaticStructuralSolver(model_part,domain_size) 
elif(Kratos_Structural_Application_var.SolverType == "DynamicSolver"):
    solver = SolverType.DynamicStructuralSolver(model_part,domain_size) 
elif(Kratos_Structural_Application_var.SolverType == "ParallelSolver"):
    solver = SolverType.StaticStructuralSolver(model_part,domain_size)
elif(Kratos_Structural_Application_var.SolverType == "ArcLengthSolver"):
    solver = SolverType.StaticStructuralSolver(model_part,domain_size)
    model_part.ProcessInfo[LAMNDA] = 0.00;  
elif(Kratos_Structural_Application_var.SolverType == "LineSearchesSolver"):
    solver = SolverType.StaticStructuralSolver(model_part,domain_size)


#solver.structure_linear_solver = Kratos_Structural_Application_var.problem_name.LinearSolver()  
if(Kratos_Structural_Application_var.LinearSolver == "SkylineLUFactorization"):
    solver.structure_linear_solver  =  SkylineLUFactorizationSolver()
elif(Kratos_Structural_Application_var.LinearSolver == "SuperLUSolver"):
    solver.structure_linear_solver  =  SuperLUSolver()
elif(Kratos_Structural_Application_var.LinearSolver == "BiCGStab_ILU0"):
    pILUPrecond = ILU0Preconditioner() 
    LST  = Kratos_Structural_Application_var.Linear_Solver_Tolerance
    LSMI = Kratos_Structural_Application_var.Linear_Solver_Max_Iteration 
    solver.structure_linear_solver  =  BICGSTABSolver(LST, LSMI, pILUPrecond)
elif(Kratos_Structural_Application_var.LinearSolver == "BiCGStab_DIAG"):
    pDiagPrecond = DiagonalPreconditioner()
    LST  = Kratos_Structural_Application_var.Linear_Solver_Tolerance
    LSMI = Kratos_Structural_Application_var.Linear_Solver_Max_Iteration  
    solver.structure_linear_solver  =  BICGSTABSolver(LST,LSMI,pDiagPrecond)
elif(Kratos_Structural_Application_var.LinearSolver == "ParallelMKLPardisoSolver"):
    solver.structure_linear_solver =  ParallelMKLPardisoSolver()

CT = Kratos_Structural_Application_var.Convergence_Tolerance;
AT = Kratos_Structural_Application_var.Absolute_Tolerance; 

if(Kratos_Structural_Application_var.Convergence_Criteria == "Displacement_Criteria"):
    solver.conv_criteria  =  DisplacementCriteria(CT,AT)
elif(Kratos_Structural_Application_var.Convergence_Criteria == "Residual_Criteria"): 
    solver.conv_criteria  =   ResidualCriteria(CT,AT)
elif(Kratos_Structural_Application_var.Convergence_Criteria == "And_Criteria"): 
    Displacement   =   DisplacementCriteria(CT,AT)
    Residual       =   ResidualCriteria(CT,AT)
    solver.conv_criteria  = AndCriteria(Residual, Displacement)
elif(Kratos_Structural_Application_var.Convergence_Criteria == "Or_Criteria"): 
    Displacement   =   DisplacementCriteria(CT,AT)
    Residual       =   ResidualCriteria(CT,AT)
    solver.conv_criteria  = OrCriteria(Residual, Displacement)

import sys
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking

solver.Initialize()
(solver).SetEchoLevel(2);


print "Forcing Node List"
forcing_nodes_list = [];

for node in model_part.Nodes:
    if(node.Y > 0.249 and node.Y < 0.2501):
        forcing_nodes_list.append(node)

for node in forcing_nodes_list:
          print node.Id , " ", node.X, " ",node.Y


def IncreaseDisplacement(forcing_nodes_list,disp):
    for node in forcing_nodes_list:
        node.SetSolutionStepValue(DISPLACEMENT_Y,0,disp)
        
 
node_1 = FindNode(model_part.Nodes, 2.5, 0.25, 0.00)
solver.CalculateReactionFlag = True;

Dt         = 1;
nsteps     = 61;
Inc        = 0.01;
disp       = 0.00;

print("initializing results")

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


for inner_step in range(1,nsteps):

    time = Dt*inner_step
    model_part.CloneTimeStep(time)

    print time
    disp = -Inc*inner_step;
    print disp; 
    IncreaseDisplacement(forcing_nodes_list,disp);   
    if(inner_step > 0):
       solver.Solve()
       if (benchmarking.InBuildReferenceMode()):
	 BenchmarkCheck(time, node_1)
       else:
         BenchmarkCheck(time, node_1)
                   

    #print the results
    print " Writing Result in Gid"
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(FORCE,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)

gid_io.FinalizeResults()
print "COMPLETED ANALYSIS"      






