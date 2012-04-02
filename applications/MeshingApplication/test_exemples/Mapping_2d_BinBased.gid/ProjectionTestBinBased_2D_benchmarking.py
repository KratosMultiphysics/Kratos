#import the configuration data as read from the GiD
import pfem_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = pfem_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = pfem_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = pfem_var.kratos_path + '/applications' ##kratos_root/applications
kratos_benchmarking_path    = pfem_var.kratos_path + '/benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_PFEMApplication = True
applications_interface.Import_MeshingApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

from KratosIncompressibleFluidApplication import *
from KratosPFEMApplication import *
from KratosMeshingApplication import *
## from now on the order is not anymore crucial
##################################################################
##################################################################
import benchmarking
import math
import pfem_var_fixed

def BenchmarkCheck(time, model_part):
    dist1 = 0.0;
    dist2 = 0.0;
    dist3 = 0.0; 
    for node in model_part.Nodes:
        if(node.Id == 1243):
            dist1 = node.GetSolutionStepValue(DISTANCE)
        if(node.Id == 924):
            dist2 = node.GetSolutionStepValue(DISTANCE)
        if(node.Id == 446):
            dist3 = node.GetSolutionStepValue(DISTANCE)
    benchmarking.Output(time, "Time")
    benchmarking.Output(dist1, "distance on node #1243", 0.00001)
    benchmarking.Output(dist2, "distance on node #924", 0.00001)
    benchmarking.Output(dist3, "distance on node #446", 0.00001)

def BenchmarkCheckPfem(time, model_part):
    dist1 = 0.0;
    dist2 = 0.0;
    dist3 = 0.0; 
    for node in model_part.Nodes:
        if(node.Id == 122):
            dist1 = node.GetSolutionStepValue(TEMPERATURE)
        if(node.Id == 199):
            dist2 = node.GetSolutionStepValue(TEMPERATURE)
        if(node.Id == 61):
            dist3 = node.GetSolutionStepValue(TEMPERATURE)
    benchmarking.Output(time, "Time")
    benchmarking.Output(dist1, "distance on node #122", 0.00001)
    benchmarking.Output(dist2, "distance on node #199", 0.00001)
    benchmarking.Output(dist3, "distance on node #61", 0.00001)
    
#defining a model part
pfem_model_part = ModelPart("FluidPart");  
fixed_model_part = ModelPart("FixedFluidPart");

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import pfem_solver_ale
pfem_solver_ale.AddVariables(pfem_model_part)
pfem_model_part.AddNodalSolutionStepVariable(DISTANCE)
pfem_model_part.AddNodalSolutionStepVariable(TEMPERATURE)

pfem_solver_ale.AddVariables(fixed_model_part)
fixed_model_part.AddNodalSolutionStepVariable(DISTANCE)
fixed_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
#reading a model
name_pfem  = pfem_var.problem_name
name_fixed = pfem_var_fixed.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

#reading the pfem model part
gid_io = GidIO(name_pfem,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_pfem = ModelPartIO(name_pfem)
model_part_io_pfem.ReadModelPart(pfem_model_part)
#reading the fixed model part
model_part_io_fixed = ModelPartIO(name_fixed)
model_part_io_fixed.ReadModelPart(fixed_model_part)
#the buffer size should be set up here after the mesh is read for the first time
pfem_model_part.SetBufferSize(2)

#importing the solver files
pfem_solver_ale.AddDofs(pfem_model_part)
for node in pfem_model_part.Nodes:
    node.AddDof(DISTANCE);
    node.AddDof(TEMPERATURE)
pfem_solver_ale.AddDofs(fixed_model_part)
for node in fixed_model_part.Nodes:
    node.AddDof(DISTANCE);
    node.AddDof(TEMPERATURE)
##print "LINE 80"
#setting the limits of the bounding box
box_corner1 = Vector(3); box_corner1[0]=-3.0; box_corner1[1]=-3.0; box_corner1[2]=-3.0;
box_corner2 = Vector(3); box_corner2[0]=20.0; box_corner2[1]=20.0;  box_corner2[2]=20.0;


Projection = BinBasedMeshTransfer2D()
node_locator = BinBasedFastPointLocator2D(fixed_model_part)
node_locator_from_elem = BinBasedNodesInElementLocator2D(fixed_model_part)

hmin = 0.4

node_locator.UpdateSearchDatabaseAssignedSize(hmin)
node_locator_from_elem.UpdateSearchDatabaseAssignedSize(hmin)

R =1.0
   
Dt = 0.1
nsteps = 100
output_Dt = 0.1
output = 0.0
output_Dt_old=0.0

tmax = 10.00

time = 0.0
step = 0
#for step in range(0,nsteps):
while (time < tmax):
    xc = time; #valid only because v ==1 m/s
    yc = 0.0
    zc = 0.0

    
    for node in pfem_model_part.Nodes:
        r2 = (node.X -xc)*(node.X -xc)+(node.Y -yc)*(node.Y -yc)
        r = math.sqrt(r2)
        if (r >= R):
            node.SetSolutionStepValue(DISTANCE,0,0.0)
        elif(r < R):
            node.SetSolutionStepValue(DISTANCE,0,R-r)

    Projection.MappingFromMovingMesh_VariableMeshes_ScalarVar(pfem_model_part,fixed_model_part, DISTANCE, DISTANCE, node_locator_from_elem);
    BenchmarkCheck(time, fixed_model_part)  

    for node in fixed_model_part.Nodes:
        node.SetSolutionStepValue(TEMPERATURE,0, node.X + time)

    Projection.DirectScalarVarInterpolation(fixed_model_part,pfem_model_part, TEMPERATURE, TEMPERATURE, node_locator)
    BenchmarkCheckPfem(time, pfem_model_part)  

    output = time - output_Dt_old

    if(output >= output_Dt):
        ##a change in the output name is needed!!!!
        res_name1 = str(name_fixed)
        gid_io.ChangeOutputName(res_name1)
        gid_io.InitializeMesh( time );
        gid_io.WriteMesh(fixed_model_part.GetMesh() )
        gid_io.FinalizeMesh();

        gid_io.InitializeResults( time, fixed_model_part.GetMesh() )
        gid_io.WriteNodalResults(DISTANCE,fixed_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TEMPERATURE,fixed_model_part.Nodes,time,0)

        gid_io.Flush()
        gid_io.FinalizeResults()

        res_name2 = str(name_pfem)
        gid_io.ChangeOutputName(res_name2)
        gid_io.InitializeMesh( time );
        gid_io.WriteMesh(pfem_model_part.GetMesh() )
        gid_io.FinalizeMesh();

        gid_io.InitializeResults( time, pfem_model_part.GetMesh() )
        gid_io.WriteNodalResults(DISTANCE,pfem_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TEMPERATURE,pfem_model_part.Nodes,time,0)

        gid_io.Flush()
        gid_io.FinalizeResults()


        output_Dt_old = time
    
        print "output step finished"
    for node in pfem_model_part.Nodes:
        node.X = node.X + Dt;
        
    time = time + Dt
    step = step + 1
    
print "solution finished *************************************************"         
        

