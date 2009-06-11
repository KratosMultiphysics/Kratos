##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
import math
import matplotlib

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path

kratos_libs_path = '/home/nelson/kratos/libs' ##kratos_root/libs
kratos_applications_path = '/home/nelson/kratos/applications' ##kratos_root/applications
kratos_python_scripts_path = '/home/nelson/kratos/applications/structural_application/python_scripts'

import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_python_scripts_path)

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

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosExternalSolversApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here 
#import structural_solver_dynamic
#import structural_solver_dynamic_superlu
#import structural_solver_relaxation
##import structural_solver_static
import structural_solver_static_arc_length

structural_solver_static_arc_length.AddVariables(model_part)
#structural_solver_relaxation.AddVariables(model_part)
#structural_solver_dynamic.AddVariables(model_part)
#structural_solver_dynamic_superlu.AddVariables(model_part)

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

gid_io = GidIO("Domo",gid_mode,multifile,deformed_mesh_flag, write_conditions)
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
#structural_solver_dynamic.AddDofs(model_part)
#structural_solver_dynamic_superlu.AddDofs(model_part)
#structural_solver_relaxation.AddDofs(model_part)
structural_solver_static_arc_length.AddDofs(model_part)

#creating a fluid solver object
#solver = structural_solver_relaxation.RelaxationStructuralSolver(model_part,domain_size)
#solver = structural_solver_dynamic.DynamicStructuralSolver(model_part,domain_size)
#solver = structural_solver_dynamic_superlu.DynamicStructuralSolver(model_part,domain_size)
solver = structural_solver_static_arc_length.StaticStructuralSolver(model_part,domain_size)

#solver.structure_linear_solver =  SuperLUSolver()

## Variables Arc Length
solver.Ide                        = 125
solver.factor_delta_lmax          = 1.00
solver.max_iteration              = 500
solver.toler                      = 1.0E-10
solver.norm                       = 1.0E-7
solver.MaxLineSearchIterations    = 20
solver.tolls                      = 0.000001          
solver.amp                        = 1.618             
solver.etmxa                      = 5                 
solver.etmna                      = 0.1               
solver.CalculateReactionFlag      = True
solver.ReformDofSetAtEachStep     = True
solver.MoveMeshFlag               = True
solver.ApplyLineSearches          = True
        


solver.Initialize()
(solver).SetEchoLevel(2);

print "Forcing Node List"
forcing_nodes_list = [];

for node in model_part.Nodes:
    if(node.Y > 24.9 and node.Y < 25.01):
        forcing_nodes_list.append(node)

for node in forcing_nodes_list:
          print node.Id , " ", node.X, " ",node.Y


##model_part.Nodes[1].SetSolutionStepValue(DISPLACEMENT_Y,0,0.00);
##model_part.Nodes[1].SetSolutionStepValue(FORCE_Y,0,0.00);



def IncreasePointLoad(Load):
	#model_part.Nodes[2].SetSolutionStepValue(FORCE_Z,0,Load);
	#model_part.Nodes[3].SetSolutionStepValue(FORCE_Z,0,Load);
	model_part.Nodes[6].SetSolutionStepValue(FORCE_Z,0,Load);
	#model_part.Nodes[7].SetSolutionStepValue(FORCE_Z,0,Load);
	#model_part.Nodes[8].SetSolutionStepValue(FORCE_Z,0,Load);
	#model_part.Nodes[9].SetSolutionStepValue(FORCE_Z,0,Load);
	#model_part.Nodes[10].SetSolutionStepValue(FORCE_Z,0,Load);


def ChangeCondition(lamda):
        #new_load = model_part.Nodes[2].GetSolutionStepValue(FORCE_Z)*lamda;
	new_load_central = model_part.Nodes[6].GetSolutionStepValue(FORCE_Z)*lamda;
	#model_part.Nodes[2].SetSolutionStepValue(FORCE_Z,0,new_load)
	#model_part.Nodes[3].SetSolutionStepValue(FORCE_Z,0,new_load)
	model_part.Nodes[6].SetSolutionStepValue(FORCE_Z,0,new_load_central)
	#model_part.Nodes[7].SetSolutionStepValue(FORCE_Z,0,new_load)
	#model_part.Nodes[8].SetSolutionStepValue(FORCE_Z,0,new_load)
	#model_part.Nodes[9].SetSolutionStepValue(FORCE_Z,0,new_load)
	#model_part.Nodes[10].SetSolutionStepValue(FORCE_Z,0,new_load)


 
## DATOS DE ENTRADA        

Dt         =  1;
nsteps     =  2000;
Load_inc   = -0.25;
P_load     = Vector(nsteps);
Carga      = Vector(nsteps);
Delta_1    = Vector(nsteps);
Delta      = Vector(nsteps);
Load_central = 0.00;
Load_externa = 0.00;
model_part.ProcessInfo[LAMNDA] = 0.00;
Carga[0]   = 0.00;
P_load[0]  = 0.00;
Delta_1[0] = 0.00;
Delta[0]   = 0.00;

print("Initializing Results. Please wait.....")
gid_io.InitializeResults(mesh_name,(model_part).GetMesh())

     
for inner_step in range(1,nsteps):

    time = Dt*inner_step
    model_part.CloneTimeStep(time)   
    model_part.ProcessInfo[TIME_STEPS] = inner_step
    print time
    
    
    #Load = Calculate_Next_Step_Load(Carga,Delta_1,inner_step,Load_inc);
    Load =  Load_inc*inner_step;
    IncreasePointLoad(Load);
    solver.Solve()

    print model_part.ProcessInfo[LAMNDA]
    Carga[inner_step]    =   -model_part.Nodes[6].GetSolutionStepValue(FORCE_Z)*model_part.ProcessInfo[LAMNDA];
    Delta[inner_step]    =   -(model_part.Nodes[6].GetSolutionStepValue(DISPLACEMENT_Z));
    P_load[inner_step]   =   -model_part.Nodes[10].GetSolutionStepValue(FORCE_Z)*model_part.ProcessInfo[LAMNDA];
    Delta_1[inner_step]  =   -(model_part.Nodes[10].GetSolutionStepValue(DISPLACEMENT_Z));
    ChangeCondition(model_part.ProcessInfo[LAMNDA])  
           
    #print the results
    print " Writing Result in Gid"
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(FORCE,model_part.Nodes,time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
    gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
    gid_io.PrintOnGaussPoints(DAMAGE,model_part,time)


print P_load
print Delta_1

from pylab import *
font = {'fontname'   : 'Courier',
        'color'      : 'r',
        'fontweight' : 'bold',
        'fontsize'   :  11}
title('SERIE C', font, fontsize=12)
#subplot(211)
#plot(Delta_1,P_load,'g-o')
#text(0.5, 2.5, 'a line', font, color='k')
#xlim(-5,5)
#ylim(-1000,1000)
#grid(True)
#xlabel('Displacement (mm)', font)
#ylabel('P Load (KN)', font)

#subplot(212)
plot(Delta,Carga,'g-o')
text(0.5, 2.5, 'a line', font, color='k')
xlim(0,12)
ylim(-800,800)
grid(True)
xlabel('Displacement (mm)', font)
ylabel('P Load (KN)', font)
#show()
savefig('FIGURE 3')

gid_io.FinalizeResults()
print "COMPLETED ANALYSIS"   
