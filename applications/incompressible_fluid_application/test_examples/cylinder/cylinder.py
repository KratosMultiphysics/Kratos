##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

print "aaa"
#importing Kratos main library
from Kratos import *
print "bbbb"
kernel = Kernel()   #defining kernel
#print kernel
print "ccc"
#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_MeshingApplication = False
applications_interface.Import_ULFApplication = False
applications_interface.ImportApplications(kernel, kratos_applications_path)
print kernel
## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosIncompressibleFluidApplication import *
##from KratosR1MeshingApplication import *

#defining a model part
print "before creation of the model part"
model_part = ModelPart("FluidPart");
print "after creation of the model part"

##importing the solver files and adding the variables
import monolithic_solver_ALE
monolithic_solver_ALE.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(ADVPROJ)
model_part.AddNodalSolutionStepVariable(DIVPROJ)
model_part.AddNodalSolutionStepVariable(THAWONE)
model_part.AddNodalSolutionStepVariable(THAWTWO)

#adding of Variables to Model Part should be here when the "very fix container will be ready"

#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("cylinder",gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(model_part)

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

##add Degrees of Freedom to all of the nodes
monolithic_solver_ALE.AddDofs(model_part)



#creating a fluid solver object
fluid_solver = monolithic_solver_ALE.MonolithicSolver(model_part,domain_size)
#pILUPrecond = ILU0Preconditioner() 
#fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
##pDiagPrecond = DiagonalPreconditioner()
##fluid_solver.velocity_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
##fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();

##for node in model_part.Nodes:
####    node.SetSolutionStepValue(IS_FLUID,0,1.0)
##    node.SetSolutionStepValue(REACTION_X,0,0.0)

fluid_solver.max_iter = 5
fluid_solver.conv_criteria = UPCriteria(1e-2,1e-4,1e-2,1e-4)
fluid_solver.Initialize()

#settings to be changed
Re = 100
nsteps =1200
output_step = 1

Dt = 1000.0/Re
if(Dt > 0.000000001):
    Dt =1
out = 0
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re)
    node.SetSolutionStepValue(DENSITY,0,1.0)
    node.SetSolutionStepValue(BODY_FORCE_X,0,0)
    node.SetSolutionStepValue(BODY_FORCE_Y,0,0)

##for node in model_part.Nodes:
##    if(node.X < 0.001 and node.Y<0.001):
##        node.Fix(PRESSURE)


zero = Vector(3);
zero[0] = 0.0;
zero[1] = 0.0;
zero[2] = 0.0;

#gid_io.InitializeResults( 0.0, model_part.GetMesh() )
firstT = 0
df = 0

for step in range(1,nsteps):
    print "line49"

    if(step>100):
        firstT=100*1
        df = 100
        Dt=.1   
        fluid_solver.max_iter =5
        fluid_solver.conv_criteria = UPCriteria(1e-6,1e-9,1e-6,1e-6)
        fluid_solver.Initialize()
        
        
    time = firstT+Dt*(step-df)
    print time
    model_part.CloneTimeStep(time)

    print "qui"

    print time
    
    #solving the fluid problem
    if(step > 3):
        fluid_solver.Solve()
        print "After solve"
        
    print "li"


    #print the results
    if(out == output_step):
        file_name = "twofluid"
        file_name = file_name + str(step)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((model_part).GetMesh());
        gid_io.WriteMesh((model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MESH_VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_STRUCTURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_BOUNDARY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VISCOSITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DIVPROJ,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(ADVPROJ,model_part.Nodes,time,0)
       # gid_io.PrintOnGaussPoints(THAWONE,model_part,time)
        #gid_io.PrintOnGaussPoints(THAWTWO,model_part,time)
        gid_io.FinalizeResults();
        
        out = 0
    out = out + 1


          
        

