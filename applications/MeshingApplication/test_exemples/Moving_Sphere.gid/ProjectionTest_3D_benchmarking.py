##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs/' ##kratos_root/libs
kratos_applications_path = '../../../../applications/' ##kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
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

from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
from KratosMeshingApplication import *
## from now on the order is not anymore crucial
##################################################################
##################################################################
import benchmarking
import math

def BenchmarkCheck(time, model_part):
    dist1 = 0.0;
    dist2 = 0.0;
    dist3 = 0.0;
    for node in model_part.Nodes:
        if(node.Id == 11241):
            dist1 = node.GetSolutionStepValue(DISTANCE)
        if(node.Id == 7742):
            dist2 = node.GetSolutionStepValue(DISTANCE)
        if(node.Id == 4535):
            dist3 = node.GetSolutionStepValue(DISTANCE)
    benchmarking.Output(time, "Time")
    benchmarking.Output(dist1, "distance on node #11241", 0.00001)
    benchmarking.Output(dist2, "distance on node #7742", 0.00001)
    benchmarking.Output(dist3, "distance on node #4535", 0.00001)

#defining a model part
pfem_model_part = ModelPart("FluidPart");  
fixed_model_part = ModelPart("FixedFluidPart");

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import pfem_solver_ale
pfem_solver_ale.AddVariables(pfem_model_part)
pfem_model_part.AddNodalSolutionStepVariable(DISTANCE)
pfem_solver_ale.AddVariables(fixed_model_part)
fixed_model_part.AddNodalSolutionStepVariable(DISTANCE)
#reading a model
name_pfem = "Moving_Sphere"
name_fixed = "Sphere_fixed"

gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

#reading the fixed model part
gid_io = GidIO(name_pfem,gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(pfem_model_part)

#reading the fixed model part
data_io = DatafileIO(name_fixed)
data_io.ReadModelPart(fixed_model_part)

#the buffer size should be set up here after the mesh is read for the first time
pfem_model_part.SetBufferSize(2)

#importing the solver files
pfem_solver_ale.AddDofs(pfem_model_part)
for node in pfem_model_part.Nodes:
    node.AddDof(DISTANCE);
    
#for node in pfem_model_part.Nodes:
#    node.GetSolutionStepValue(TEMPERATURE);

#setting the limits of the bounding box
box_corner1 = Vector(3); box_corner1[0]=-3.0; box_corner1[1]=-3.0; box_corner1[2]=-3.0;
box_corner2 = Vector(3); box_corner2[0]=20.0; box_corner2[1]=20.0;  box_corner2[2]=20.0;

#creating a fluid solver object
name = str("Sphere_pfem")


##for node in pfem_model_part.Nodes:
##    node.SetSolutionStepValue(VELOCITY_X,0,1.0)
##    node.SetSolutionStepValue(VISCOSITY,0,0.000001)
##    node.SetSolutionStepValue(DENSITY,0,1000.000)

Projection = MeshTransfer3D()

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
        r2 = (node.X -xc)*(node.X -xc)+(node.Y -yc)*(node.Y -yc)+(node.Z -zc)*(node.Z -zc)
        r = math.sqrt(r2)
        if (r >= R):
            node.SetSolutionStepValue(DISTANCE,0,0.0)
        elif(r < R):
            node.SetSolutionStepValue(DISTANCE,0,R-r)

    Projection.DirectScalarVarInterpolation(pfem_model_part,fixed_model_part, DISTANCE);
    BenchmarkCheck(time, fixed_model_part)  

    output = time - output_Dt_old
##    if(step == 0):
##        res_name1 = str("Mesh_fixed")
##        gid_io.ChangeOutputName(res_name1)
##        gid_io.InitializeMesh( time );
##        gid_io.WriteMesh(fixed_model_part.GetMesh() )
##        gid_io.FinalizeMesh();
##
##        gid_io.InitializeResults( time, fixed_model_part.GetMesh() )
##        gid_io.WriteNodalResults(DISTANCE,fixed_model_part.Nodes,time,0)
##        gid_io.FinalizeResults()
##
##        res_name2 = str("Mesh_pfem")
##        gid_io.ChangeOutputName(res_name2)
##        gid_io.InitializeMesh( time );
##        gid_io.WriteMesh(pfem_model_part.GetMesh() )
##        gid_io.FinalizeMesh();
##
##        gid_io.InitializeResults( time, fixed_model_part.GetMesh() )
##        gid_io.WriteNodalResults(DISTANCE,pfem_model_part.Nodes,time,0)
##        gid_io.FinalizeResults()
##
##        gid_io.Flush()


    if(output >= output_Dt):
        ##a change in the output name is needed!!!!
        res_name1 = str("Mesh_fixed")
        gid_io.ChangeOutputName(res_name1)
        gid_io.InitializeMesh( time );
        gid_io.WriteMesh(fixed_model_part.GetMesh() )
        gid_io.FinalizeMesh();

        gid_io.InitializeResults( time, fixed_model_part.GetMesh() )
        gid_io.WriteNodalResults(DISTANCE,fixed_model_part.Nodes,time,0)

        gid_io.Flush()
        gid_io.FinalizeResults()

        res_name2 = str("Mesh_pfem")
        gid_io.ChangeOutputName(res_name2)
        gid_io.InitializeMesh( time );
        gid_io.WriteMesh(pfem_model_part.GetMesh() )
        gid_io.FinalizeMesh();

        gid_io.InitializeResults( time, fixed_model_part.GetMesh() )
        gid_io.WriteNodalResults(DISTANCE,pfem_model_part.Nodes,time,0)

        gid_io.Flush()
        gid_io.FinalizeResults()


        output_Dt_old = time
    
        print "output step finished"
    for node in pfem_model_part.Nodes:
        node.X = node.X + Dt;
        
    time = time + Dt
    step = step + 1
    
print "solution finished *************************************************"         
        

