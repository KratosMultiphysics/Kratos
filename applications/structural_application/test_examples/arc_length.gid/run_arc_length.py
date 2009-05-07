##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path

kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
kratos_python_scripts_path = '../../../../applications/structural_application/python_scripts'
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking


import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_python_scripts_path)
sys.path.append(kratos_benchmarking_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_StructuralApplication = True
applications_interface.Import_KratosExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *
print kernel


import benchmarking

def FindNode(node_list,x,y,z):
    for node in node_list:
        if ((node.X - x) ** 2 + (node.Y - y) ** 2 + (node.Z - z) ** 2 < 0.0000001):
            print  node
            return node
    
def BenchmarkCheck(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_Y), "Node 1 Desplacement_y", 1.0)
    benchmarking.Output(node1.GetSolutionStepValue(FORCE_Y), "Node 1 Force_y", 1.0)

def AnalyticalResults(time, node1):
    benchmarking.Output(time, "Time")
    benchmarking.Output(-5.07526446248e-06 , "Node 1 Displacement_y", 1.0)
    benchmarking.Output(-0.999984773962 , "Node 1 Force_y", 1.0)


## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosExternalSolversApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#import structural_solver_static_arc_length
import structural_solver_static

#structural_solver_static_arc_length.AddVariables(model_part)
structural_solver_static.AddVariables(model_part)

###################################################################
###################################################################

model_part.AddNodalSolutionStepVariable(FORCE);
model_part.AddNodalSolutionStepVariable(DENSITY);
model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
model_part.AddNodalSolutionStepVariable(REACTION);

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
#gid_mode = GiDPostMode.GiD_PostAscii
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

gid_io = GidIO("arc_length",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print model_part

#writing the mesh
#gid_io.WriteMesh(model_part.GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
#structural_solver_static_arc_length.AddDofs(model_part)
structural_solver_static.AddDofs(model_part)

#creating a fluid solver object
#solver = structural_solver_static_arc_length.StaticStructuralSolver(model_part,domain_size)
solver = structural_solver_static.StaticStructuralSolver(model_part,domain_size)

node_1 = FindNode(model_part.Nodes, 2.5, 0.25, 0.00)

## Variables Arc Length
solver.Ide               = 5
solver.factor_delta_lmax = 1.00
solver.max_iteration     = 20
solver.toler             = 1.0E-10
solver.norm              = 1.0E-7


solver.Initialize()
(solver).SetEchoLevel(2);

print "Forcing Node List"
forcing_nodes_list = [];

for node in model_part.Nodes:
    if(node.Y > 0.249 and node.Y < 0.2501):
        forcing_nodes_list.append(node)

for node in forcing_nodes_list:
          print node.Id , " ", node.X, " ",node.Y


model_part.Nodes[1].SetSolutionStepValue(FORCE_Y,0,0.00);



def IncreasePointLoad(forcing_nodes_list,Load):
    for node in forcing_nodes_list:
        node.SetSolutionStepValue(FORCE_Y,0,Load)



def Signo(m):
    if (m<0.00):
      return  -1.00;
    elif (m>0.00):
      return 1.00;
    else:
      return 0.00;


def ChangeCondition(forcingnodelist,lamda):
    for node in forcing_nodes_list:
        new_load = node.GetSolutionStepValue(FORCE_Y)*lamda;
        node.SetSolutionStepValue(FORCE_Y,0,new_load)
 

## DATOS DE ENTRADA        

Dt         =  1;
nsteps     =  2;
Load_inc   = -1;
P_load     = Vector(nsteps);
Carga      = Vector(nsteps);
Delta_1    = Vector(nsteps);
Load       = 0.00;
model_part.ProcessInfo[LAMNDA] = 0.00;
Carga[0]   = 0.00;
P_load[0]  = 0.00;
Delta_1[0] = 0.00;


print("Initializing Results. Please wait.....")
gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

     
for inner_step in range(1,nsteps):

    time = Dt*inner_step
    model_part.CloneTimeStep(time)   
    model_part.ProcessInfo[TIME_STEPS] = inner_step
    print time
    
    
    Load =  -0.999984773962 
    IncreasePointLoad(forcing_nodes_list,Load);
    
    if(inner_step > 0):
       solver.Solve()
       if (benchmarking.InBuildReferenceMode()):
	 BenchmarkCheck(time, node_1)
       else:     
         AnalyticalResults(time, node_1)
           
    #print the results
    print " Writing Result in Gid"
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(FORCE,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
    gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
    gid_io.PrintOnGaussPoints(DAMAGE,model_part,time)


gid_io.FinalizeResults()
print "COMPLETED ANALYSIS"   
