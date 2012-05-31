#import the configuration data as read from the GiD
import DEM_FEM_Application_var

##################################################################
from KratosMultiphysics import *
from KratosMultiphysics.DEM_FEM_Application import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.DEMApplication import *
##################################################################


#defining a model part
model_part = ModelPart("DEM_FEM_Part");  
model_part.AddNodalSolutionStepVariable(FORCE);


#adding of Variables to Model Part should be here when the "very fix container will be ready"
import DEM_FEM_Explicit_Solve_Strategy as DEM_FEM_Explicit_Solve_Strategy
DEM_FEM_Explicit_Solve_Strategy.AddVariables(model_part)
 

#reading a model
domain_size = DEM_FEM_Application_var.domain_size
name = DEM_FEM_Application_var.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteSphereMesh(model_part.GetMesh() )
gid_io.WriteMesh((model_part).GetMesh());

gid_io.FinalizeMesh()

print model_part
print model_part.Properties


#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)


#importing the solver files
DEM_FEM_Explicit_Solve_Strategy.AddDofs(model_part)
solver = DEM_FEM_Explicit_Solve_Strategy.DynamicStructuralSolver(model_part,domain_size)


##choosing the default value for the constitutive law 
if(DEM_FEM_Application_var.ConstitutiveLaw == "LinearElastic"):
    if(domain_size == 2):
        for prop in model_part.Properties:
            prop.SetValue(CONSTITUTIVE_LAW, Isotropic2D() )
    else:
        for prop in model_part.Properties:
            prop.SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
elif(DEM_FEM_Application_var.ConstitutiveLaw == "DruckerPrager"):    
    if(domain_size == 3):
        for prop in model_part.Properties:
            prop.SetValue(CONSTITUTIVE_LAW, DruckerPrager())
    else:
        print "Error!!!!!DruckerPrager Constitutive Law only suitable for 3D Mesh"


####choosing the liear_solver
solver.structure_linear_solver  =  SkylineLUFactorizationSolver()


####initial the solver parameter

solver.damping_ratio            = DEM_FEM_Application_var.DampRatio
solver.max_delta_time           = DEM_FEM_Application_var.Time_Step
solver.ConvUnbalForceRatio      = DEM_FEM_Application_var.ConvUnBalForceRatio
solver.contact_stiffness_ratio  = DEM_FEM_Application_var.ContactStiffRatio
solver.search_radius_extension  = DEM_FEM_Application_var.ParticleSearchTolerance

gravity = Vector(3)
gravity[0] = DEM_FEM_Application_var.gravity_x 
gravity[1] = DEM_FEM_Application_var.gravity_y 
gravity[2] = DEM_FEM_Application_var.gravity_z 

solver.gravity = gravity

if(DEM_FEM_Application_var.MassType == "VirtualMass"):
    solver.virtual_mass = True
else:
    solver.virtual_mass = False

if(DEM_FEM_Application_var.DampType == "LocalDamp"):
    solver.damp_type = 1
else:
    solver.damp_type = 2

if(DEM_FEM_Application_var.ComputeFemFemContact == "True"):
    solver.ComputeFemFemContact = True
else:
    solver.ComputeFemFemContact = False

if(DEM_FEM_Application_var.ComputeParticleRotation == "True"):
    solver.Particle_If_Cal_Rotate = 1
else:
    solver.Particle_If_Cal_Rotate = 0

if(DEM_FEM_Application_var.ComputeParticleRotationSpring == "True"):
    solver.Particle_If_Cal_Rotate_Spring = 1
else:
    solver.Particle_If_Cal_Rotate_Spring = 0

if(DEM_FEM_Application_var.ParticleNormalForceCalType == "Incremental"):
    solver.force_calculation_type_id = 1
else:
    solver.force_calculation_type_id = 2

if(DEM_FEM_Application_var.IfContinuousProblem == "True"):
    solver.continuum_simulating_OPTION = True
else:
    solver.continuum_simulating_OPTION = False

if(DEM_FEM_Application_var.IfCalInitialDelt == "True"):
    solver.delta_OPTION = True
else:
    solver.delta_OPTION = False

## Creat the expicit object

solver.Initialize()

(solver.solver).SetEchoLevel(2);

gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


Nsteps         = DEM_FEM_Application_var.TotalSteps
OutputInterval =  DEM_FEM_Application_var.OutputInterval


for step in range(1,Nsteps):
    model_part.ProcessInfo[TIME_STEPS] = step
    if(step > 0):       
        solver.Solve()

        time = model_part.ProcessInfo[TIME]
        if(step % OutputInterval == 0):
            gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)
	    ##gid_io.WriteNodalResults(FORCE,model_part.Nodes,time,0)
	    gid_io.WriteNodalResults(REACTION,model_part.Nodes,time,0)
	    gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
	    gid_io.WriteNodalResults(RHS,model_part.Nodes,time,0)
            gid_io.WriteNodalResults(RADIUS, model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(PARTICLE_NUMBER_OF_NEIGHBOURS, model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(PARTICLE_FAILURE_ID, model_part.Nodes, time, 0)
            gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR,model_part,time)
            gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR,model_part,time)
            if(solver.Particle_If_Cal_Rotate == 1):
                gid_io.WriteNodalResults(PARTICLE_ROTATION_ANGLE, model_part.Nodes, time, 0)
             


        unbalratio = model_part.ProcessInfo[DEM_FEM_CONVERGENCE_RATIO]
        if(unbalratio < solver.ConvUnbalForceRatio):            
            print "***************************Congratulations!**************************"
            print "**********Numerical System has reached a convergence state!**********"
            print "**************the loop has been stoped ahead of time!!***************"
            break
        elif (unbalratio > 1.0e10 or unbalratio < -1.0e10):
            print "***************ERROR!!!!!:The Numerical System is Out of Converginece**************"
            print "***************Please Check the Time Step or Material of the Model*****************"
            break
                     

print "Analysis Completed "

gid_io.FinalizeResults()
          
        

