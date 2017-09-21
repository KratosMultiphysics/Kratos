import fluid_only_var
import fluid_ulf_var
##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

#importing applications
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");
reduced_model_part = ModelPart("ReducedModelPart") 

aux_model_part = ModelPart("AuxModelPart")
pseudo_lag_part = ModelPart("PseudoLagPart")

Qcomp_model_part = ModelPart("pFluidPart");
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");

particle_model_part = ModelPart("FixedFluidPart")

cut_model_part = ModelPart("CutPart");	

projection_conditions_model_part = ModelPart("ProjectionPart")

lagrangian_surface_model_part = ModelPart("LagrangianSurfacePart")

import ProjectParameters
##import Rad1
##import ProjectParameters_pol

import frac_step_solverQ
frac_step_solverQ.AddVariables(Qcomp_model_part,particle_model_part)

fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
fluid_model_part.AddNodalSolutionStepVariable(AUX_VEL)
fluid_model_part.AddNodalSolutionStepVariable(DISABLE)
fluid_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
fluid_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
fluid_model_part.AddNodalSolutionStepVariable(IS_FLUID)
fluid_model_part.AddNodalSolutionStepVariable(VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(NODAL_H)

Qcomp_model_part.AddNodalSolutionStepVariable(VELOCITY)
Qcomp_model_part.AddNodalSolutionStepVariable(FRACT_VEL)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_LAGRANGIAN_INLET)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
Qcomp_model_part.AddNodalSolutionStepVariable(NODAL_H)
Qcomp_model_part.AddNodalSolutionStepVariable(DISABLE)
Qcomp_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)

reduced_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)

cut_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
cut_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
################################################
################################################
pseudo_lag_part.AddNodalSolutionStepVariable(IS_INTERFACE);
pseudo_lag_part.AddNodalSolutionStepVariable(TEMPERATURE);
################################################
################################################

#introducing input file name
input_file_name = fluid_only_var.problem_name
#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
#write_conditions = WriteConditionsFlag.WriteElementsOnly
write_conditions = WriteConditionsFlag.WriteConditions 
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)



input_file_name = fluid_ulf_var.problem_name
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_origin = ModelPartIO(input_file_name)
model_part_io_origin.ReadModelPart(Qcomp_model_part)

import math

fluid_model_part.SetBufferSize(3)
reduced_model_part.SetBufferSize(3)


Qcomp_model_part.SetBufferSize(3)

##fractional_step_solver.AddDofs(fluid_model_part)

frac_step_solverQ.AddDofs(Qcomp_model_part,particle_model_part)
subdomain_disable_process=SubdomainDisableProcess()


##solver polimero
box_corner1 = Vector(3); 
box_corner1[0]= -1.5; box_corner1[1]=-0.49; box_corner1[2]=-1.5;
box_corner2 = Vector(3); 
box_corner2[0]= 1.5; box_corner2[1]=1.5; box_corner2[2]=1.5;

fluid_solverq = frac_step_solverQ.FracStepSolver(Qcomp_model_part,particle_model_part,domain_size)
fluid_solverq.laplacian_form = 1; 
fluid_solverq.predictor_corrector = True
fluid_solverq.ReformDofAtEachIteration = True
fluid_solverq.linear_solver = SkylineLUFactorizationSolver()
fluid_solverq.max_vel_its = 5
fluid_solverq.time_order = 1
fluid_solverq.echo_level = 1
fluid_solverq.Initialize()


time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
lin_solver =  SkylineLUFactorizationSolver()
ReformDofSet=True
ProjDirichletLinStrat=ResidualBasedLinearStrategy(aux_model_part, time_scheme, lin_solver, False, ReformDofSet, False, False)
ProjDirichletLinStrat.SetEchoLevel(2)

proj_dirichlet_process=ApplyProjDirichletProcess()
pseudo_lag_process=PseudoLagPartProcess

import NistParameters

neigh_finder = FindNodalNeighboursProcess(aux_model_part,9,18)
elem_neigh_f=FindElementalNeighboursProcess(aux_model_part, 2, 10)
cond_neigh_f=FindConditionsNeighboursProcess(aux_model_part, 2, 10)
AAA = MeshTransfer3D()

#mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( cut_model_part.GetMesh() ) 					#originally: gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()
gid_io.InitializeResults(mesh_name,(cut_model_part).GetMesh())			#originally: gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())
#####################################################################
EmbeddedUtils=EmbeddedUtils()
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
lin_solver =  SkylineLUFactorizationSolver()
ReformDofSet=True
ProjDirichletLinStrat=ResidualBasedLinearStrategy(projection_conditions_model_part, time_scheme, lin_solver, False, ReformDofSet, False, False)
ProjDirichletLinStrat.SetEchoLevel(2)

SaveLagrangianSurfaceProcess=SaveLagrangianSurfaceProcess()
calculate_distance_process = CalculateSignedDistanceTo3DConditionSkinProcess(lagrangian_surface_model_part,fluid_model_part);


##SaveLagrangianSurfaceProcess.SaveSurfaceConditions( Qcomp_model_part, lagrangian_surface_model_part)
##calculate_distance_process.Execute()
##
##for node in fluid_model_part.Nodes:
##    node.SetSolutionStepValue(IS_INTERFACE,0,0)
##EmbeddedUtils.CreateIntersConditions(fluid_model_part, projection_conditions_model_part)
##SubDomainDisableProcess=SubdomainDisableProcess()
#SubDomainDisableProcess.SaveReducedPart(fluid_model_part, reduced_model_part)
#if (projection_conditions_model_part.Conditions.Size()>0):
#    print ("SOLVING THE INTERFACE B.C. PROBLEM")
#    print ("SOLVING THE INTERFACE B.C. PROBLEM")
#    ProjDirichletLinStrat.Solve()
ApplyProjDirichletProcess=ApplyProjDirichletProcess()
ApplyProjDirichletProcess.ApplyProjDirichlet(fluid_model_part)
#subdomain_disable_process.SaveReducedPart(fluid_model_part, reduced_model_part)
#subdomain_disable_process.SaveReducedPart(mf_model_part, mf_model_part_red)

 
for node in Qcomp_model_part.Nodes:
    node.SetSolutionStepValue(NODAL_H,0,0.0006)  ##0.0025  ##0.0012  ##0.005 0.0009

particle_utils = ParticleUtils2D()    
           
time=0.0
final_time=10
step=0
output_step=0.001
next_output_time = 0
Dt=0.005
while(time < final_time):
         
    Dt = 0.005; 

    time = time + Dt
    fluid_model_part.CloneTimeStep(time)
    reduced_model_part.ProcessInfo=fluid_model_part.ProcessInfo
    aux_model_part.ProcessInfo=fluid_model_part.ProcessInfo
    combined_model_part.ProcessInfo=fluid_model_part.ProcessInfo
    
    Qcomp_model_part.CloneTimeStep(time)
    



    fluid_solverq.RemeshAux()	

    lagrangian_surface_model_part.Conditions.clear()
     
    SaveLagrangianSurfaceProcess.SaveSurfaceConditions(Qcomp_model_part, lagrangian_surface_model_part)


    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(IS_INTERFACE,0,0)
        node.SetSolutionStepValue(DISTANCE,0,0.0)
        node.SetSolutionStepValue(DISABLE,0,0)       
        


    calculate_distance_process = CalculateSignedDistanceTo3DConditionSkinProcess(lagrangian_surface_model_part,fluid_model_part);
    calculate_distance_process.Execute()
 
           

    
    if(time > next_output_time):
        
        res_name4 = str("POLIMERO")
        gid_io.ChangeOutputName(res_name4)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((Qcomp_model_part).GetMesh());
        gid_io.WriteMesh((Qcomp_model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(time, (Qcomp_model_part).GetMesh());
        gid_io.WriteNodalResults(IS_STRUCTURE,Qcomp_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FREE_SURFACE,Qcomp_model_part.Nodes,time,0)
        gid_io.Flush()
        gid_io.FinalizeResults();
        
        res_name4 = str("")
        gid_io.ChangeOutputName(res_name4)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((fluid_model_part).GetMesh());
        gid_io.WriteMesh((fluid_model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(time, (fluid_model_part).GetMesh());
        gid_io.WriteNodalResults(DISTANCE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_INTERFACE,fluid_model_part.Nodes,time,0)
        gid_io.Flush()
        gid_io.FinalizeResults();

        particle_utils.VisualizationModelPart(cut_model_part, fluid_model_part,Qcomp_model_part )
        
##       
        res_name4 = str("SUMA")
        gid_io.ChangeOutputName(res_name4)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((cut_model_part).GetMesh());
        gid_io.WriteMesh((cut_model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(time, (cut_model_part).GetMesh());
##        gid_io.WriteNodalResults(VELOCITY,cut_model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(YCH4,cut_model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(TEMPERATURE,cut_model_part.Nodes,time,0)
##        gid_io.WriteNodalResults(IS_FREE_SURFACE,cut_model_part.Nodes,time,0)        
        gid_io.WriteNodalResults(IS_INTERFACE,cut_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISTANCE,cut_model_part.Nodes,time,0) 
        gid_io.Flush()
        gid_io.FinalizeResults();
 
        
        next_output_time = next_output_time  + output_step;
	

        out = 0
    out = out + 1
    step = step + 1
    
