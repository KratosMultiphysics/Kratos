def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node
    
def BenchmarkCheck(time, node1, node2, node3, node4):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_X), "Node 1 Displacement_x", 0.00001)
    benchmarking.Output(node2.GetSolutionStepValue(DISPLACEMENT_Y), "Node 2 Displacement_y", 0.00001)
    benchmarking.Output(node3.GetSolutionStepValue(REACTION_X), "Node 3 Reaction_x", 0.00001)
    benchmarking.Output(node4.GetSolutionStepValue(REACTION_Y), "Node 4 Reaction_y", 0.00001)



#def AnalyticalResults(time, node1, node2,node3, node4):
    #benchmarking.Output(time, "Time")
    #benchmarking.Output(-0.221921365586,  "Node 1 Displacement_x", 0.00001)
    #benchmarking.Output(-0.0361068223759, "Node 2 Displacement_y", 0.00001)
    #benchmarking.Output( 51.6844785228,   "Node 3 Reaction_x", 0.00001)
    #benchmarking.Output( -123.134969306,   "Node 4 Reaction_y", 0.00001)

##################################################################
##################################################################

import sys
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking

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


##################################################################
##################################################################

from time import *
print ctime()
t0 = clock()

#including kratos path
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
    
    
#setting the domain size for the problem to be solved
domain_size = Kratos_Structural_Application_var.domain_size


#defining a model part
model_part = ModelPart("StructurePart");  
model_part.AddNodalSolutionStepVariable(FORCE);
if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
  model_part.AddNodalSolutionStepVariable(ROTATION);


import structural_solver_dynamic as SolverType
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

model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D())
print "Linear elastic model selected"

print model_part
print model_part.Properties

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

##importing the rotational dofs degrees of freedom if necessary
RotationalDofs()


#importing the solver files
SolverType.AddDofs(model_part)
solver = SolverType.DynamicStructuralSolver(model_part,domain_size) 
solver.structure_linear_solver  =  SuperLUSolver()
solver.CalculateReactionFlag = True;


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



solver.structure_linear_solver  =  SkylineLUFactorizationSolver()


node_1 = FindNode(model_part.Nodes, 0.05, 1.00, 0.00)
node_2 = FindNode(model_part.Nodes, 0.00, 1.00, 0.00)
node_3 = FindNode(model_part.Nodes, 0.00, 0.00, 0.00)
node_4 = FindNode(model_part.Nodes, 0.05, 0.00, 0.00)

solver.Initialize()
(solver).SetEchoLevel(2);


Dt = 0.001
nsteps =11
print("initializing results")
gid_io.InitializeResults(mesh_name,(model_part).GetMesh())
for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)
    #print model_part.ProcessInfo()[TIME]
    #solving the fluid problem
    if(step > 3):
       solver.Solve()
       if (benchmarking.InBuildReferenceMode()):
         #AnalyticalResults(time, node_1, node_2, node_3, node_4)
	 BenchmarkCheck(time, node_1, node_2, node_3, node_4)
       else:
         BenchmarkCheck(time, node_1, node_2, node_3, node_4)
    #print the results
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
gid_io.FinalizeResults()
print "Completed Analysis"


          
        