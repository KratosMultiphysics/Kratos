def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node
    
def BenchmarkCheck(time, node1, node2):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_X), "Node 1 Displacement_x", 0.0000001)
    benchmarking.Output(node2.GetSolutionStepValue(DISPLACEMENT_Y), "Node 2 Displacement_y", 0.0000001)

def AnalyticalResults(time, node1, node2):
    benchmarking.Output(time, "Time")
    benchmarking.Output(0.001000, "Node 1 Displacement_x", 0.0000001)

import sys
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking

#import the configuration data as read from the GiD
import Kratos_Structural_Application_var

from time import *
print ctime()
t0 = clock()

#including kratos path
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
    
    
#setting the domain size for the problem to be solved
domain_size = Kratos_Structural_Application_var.domain_size


#defining a model part
model_part = ModelPart("StructurePart");  
model_part.AddNodalSolutionStepVariable(FORCE);
if(Kratos_Structural_Application_var.Rotational_Dofs == "True"):
  model_part.AddNodalSolutionStepVariable(ROTATION);


import structural_solver_static as SolverType
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


model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D())
print "Linear elastic model selected"

print model_part
print model_part.Properties

#writing the mesh
#gid_io.WriteMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
SolverType.AddDofs(model_part)


#creating a fluid solver object
#solver = structural_solver_relaxation.RelaxationStructuralSolver(model_part,domain_size)
#solver = structural_solver_dynamic.DynamicStructuralSolver(model_part,domain_size)
#solver = structural_solver_dynamic_superlu.DynamicStructuralSolver(model_part,domain_size)
solver = SolverType.StaticStructuralSolver(model_part,domain_size)
node_1 = FindNode(model_part.Nodes, 1.00, 0.00, 0.00)
node_2 = FindNode(model_part.Nodes, 0.165, 0.745, 0.702)


solver.Initialize()
(solver).SetEchoLevel(2);


Dt         = 1;
nsteps     = 4;

print("initializing results")

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    if(step>2):
        solver.Solve()
	
	if (benchmarking.InBuildReferenceMode()):

	  AnalyticalResults(time, node_1, node_2)
	else:

	  BenchmarkCheck(time, node_1, node_2)

print "COMPLETED ANALYSIS"   
