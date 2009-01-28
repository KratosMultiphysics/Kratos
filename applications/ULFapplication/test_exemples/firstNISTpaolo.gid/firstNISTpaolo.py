#domain size
domain_size = 2
#total simulation time
total_time = 2.0
#the max time step - it may be decreased in case it is necessary to avoid element inversion
max_delta_time = 0.001
#output time (every xxx seconds)
output_dt = 0.005
#safety factor for the delta time estimation
safety_factor = 0.5
#PATH of where the kratos library is installed
kratos_libs_path = '../../../../libs/'
#PATH of where your application is installed
kratos_applications_path = '../../../../applications/'
project_name = 'firstNISTpaolo'
##################################################################
##################################################################
## ATTENTION: here the order is important

import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel


#importing applications
import applications_interface
applications_interface.Import_PFEMApplication = True
applications_interface.Import_ULFApplication = True
applications_interface.Import_StructuralApplication = True
applications_interface.Import_ConvectionDiffusionApplication = True
###applications_interface.Import_ExternalSolversApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)


## from now on the order is not anymore crucial
##################################################################
##################################################################
print kernel
print "aaa"
from KratosULFApplication import *
print "bbb"
##from KratosPFEMApplication import *
from KratosConvectionDiffusionApplication import *
print "ccc"
from KratosExternalSolversApplication import *

def PrintLevel(time,outfile,level):
    out = str(time) + " " + str(level) + "\n"
    outfile.write(out)
    outfile.flush()

#Calculate mass loss rate -- assume sample is 10 cm wide
def PrintVolume(time,outfile,volume0):
    checkresults = solver.CheckForInvertedElements()
    volume = checkresults.pop()
    volume_percent = (volume/volume0)*100.
    density = 900. # (kg/m^3) Should actually use value set in NISTParameters.py
    thickfrac = 0.1 # fraction of 1 meter for thickness perpendicular to 2D plane
    mass = volume * thickfrac * density * 1000. # convert to grams
    outstring = str(time) + " "
    outstring += str(volume) + " "
    outstring += str(volume_percent) + " "
    outstring += str(mass) + "\n"
    outfile.write( outstring )
    outfile.flush()

#Calculate mass above and below y=0 line -- assume sample is 10 cm wide
def PrintVolumeCatchpan(time,outfile,volume0):
    checkresults = solver.CheckForInvertedElements()
    volume = checkresults.pop()
    volume_percent = (volume/volume0)*100.
    objectvol = 0.
    catchpanvol = 0.
    for node in fluid_model_part.Nodes:
        if (node.Y > 0.0):
            objectvol=objectvol + node.GetSolutionStepValue(NODAL_AREA)
        else:
            catchpanvol=catchpanvol + node.GetSolutionStepValue(NODAL_AREA)
            
    density = 900. # (kg/m^3) Should actually use value set in NISTParameters.py
    thickfrac = 0.1 # fraction of 1 meter for thickness perpendicular to 2D plane
    mass = volume * thickfrac * density * 1000. # convert to grams
    objectmass = objectvol * thickfrac * density * 1000. # convert to grams
    catchpanmass = catchpanvol * thickfrac * density * 1000. # convert to grams
    outstring = str(time) + " "
    outstring += str(objectvol) + " "
    outstring += str(objectmass) + " "
    outstring += str(catchpanvol) + " "
    outstring += str(catchpanmass) + " "
    outstring += str(volume) + " "
    outstring += str(mass) + " "
    outstring += str(volume_percent) + "\n"
    outfile.write( outstring )
    outfile.flush()

#this is for the debugging
##x = raw_input("stopped to allow debug: set breakpoints and press enter to continue");

#defining a model part
fluid_model_part = ModelPart("FluidPart");  
temperature_model_part = ModelPart("TemperaturePart");  
structure_model_part = ModelPart("StructurePart");  
combined_model_part = ModelPart("CombinedPart");

#adding of Variables to Model Part should be here when the "very fix container will be ready"
#importing the solver files
import ulf_fsi
ulf_fsi.AddVariables(fluid_model_part)

import nonlinear_convection_diffusion_solver
nonlinear_convection_diffusion_solver.AddVariables(fluid_model_part) #the nodes are the same 

#reading a model
gid_io = GidIO(project_name,GiDPostMode.GiD_PostBinary)

gid_io.ReadModelPart(fluid_model_part)
#gid_io.ReadMesh(model_part.GetMesh())
print fluid_model_part

for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(BULK_MODULUS,0,-10000.0);

fluid_model_part.GetMesh().Properties[1][BULK_MODULUS] = -10000.0


#generating temperature model part
from KratosIncompressibleFluidApplication import *
##from KratosPFEMApplication import *
from KratosULFApplication import *
#import KratosPFEMApplication
NISTTools = NistUtils()
NISTTools.GenerateModelPart(fluid_model_part,temperature_model_part,domain_size);

#the buffer size should be set up here after the mesh is read for the first time
fluid_model_part.SetBufferSize(2)

ulf_fsi.AddDofs(fluid_model_part)

nonlinear_convection_diffusion_solver.AddDofs(temperature_model_part) #the nodes are the same 

import NistParameters
NistParameters.InitialConditions(fluid_model_part)

#setting the limits of the bounding box
box_corner1 = Vector(3); box_corner1[0]=-0.1030001; box_corner1[1]=-0.0330001; box_corner1[2]=-0.1;
box_corner2 = Vector(3); box_corner2[0]=0.1030001; box_corner2[1]=0.100001;  box_corner2[2]=0.1;



#creating a fluid solver object
name = project_name
solver = ulf_fsi.ULF_FSISolver(fluid_model_part, structure_model_part, combined_model_part, box_corner1, box_corner2, domain_size)
solver.alpha_shape = 1.5;
solver.echo_level = 1;
solver.h_multiplier = 0.3
solver.model_linear_solver = SkylineLUFactorizationSolver()  # not available on Windows
#solver.model_linear_solver = SuperLUSolver()  # not available on Windows
#initializing the solver
solver.Initialize()



scalarout = open("volume_history.out", 'w')
scalarout.write("time  volume(m^2)  volumepercent(%)  mass(g) \n")
checkresults = solver.CheckForInvertedElements()
volume0 = checkresults.pop()

catchpanout = open("volume_history_catchpan.out", 'w')
catchpanout.write("time(s) objectvol(m^2) objectmass(g) catchpanvol(m2) catchpanmass(g) totalvol(m2) totalmass(g) totalvolpercent(%) \n")
checkresults = solver.CheckForInvertedElements()
volume0 = checkresults.pop()

#convection diffusion solver
temperature_solver = nonlinear_convection_diffusion_solver.ConvectionDiffusionSolver(temperature_model_part,domain_size)
#temperature_solver = convection_diffusion_solver.ConvectionDiffusionSolver(temperature_model_part,domain_size)
temperature_solver.time_order = 1
temperature_solver.ReformDofAtEachIteration = True
temperature_solver.echo_level = 0
temperature_solver.Initialize()


Dt = 0.001
nsteps = 20000
#output_Dt = 0.05
output_Dt = 0.5
min_dt = 0.0002
#max_dt = 0.05
max_dt = 0.05
safety_factor = 0.1

next_output_time = output_Dt

time = Dt
#OutputStep(time,gid_io,fluid_model_part,domain_size)


NistParameters.InitialConditions(fluid_model_part)

NISTTools.ApplyInitialTemperature(fluid_model_part,500.0);
NistParameters.CalculateViscosity(fluid_model_part.Nodes)


face_heat_util = FaceHeatUtilities()

for node in fluid_model_part.Nodes:
        node.Fix(DISPLACEMENT_Z);
        node.SetSolutionStepValue(DISPLACEMENT_Z,0,0.0);

time = 0.00
step = 0
PrintVolume (time, scalarout, volume0)
PrintVolumeCatchpan (time, catchpanout, volume0)

while(time < 1500.0):
    step = step + 1

    print "min_dt", min_dt
    print "max_dt", max_dt
    new_Dt = solver.EstimateDeltaTime(max_dt,domain_size)
    print "forever dt", new_Dt
    if(step < 10):
        new_Dt = 0.01*new_Dt


    if(time < 150):
        new_Dt =4.0



    time = time + new_Dt*safety_factor

    fluid_model_part.CloneTimeStep(time)
    structure_model_part.CloneTimeStep(time)
    combined_model_part.CloneTimeStep(time)
    temperature_model_part.CloneTimeStep(time)

    print time

    #solving the fluid problem
    if(step > 3):

##        NistParameters.ApplyBoundaryconditions(fluid_model_part.Nodes)
        NISTTools.GenerateModelPart(fluid_model_part,temperature_model_part,domain_size);
        face_heat_util.ApplyFaceHeat(fluid_model_part.Conditions,30000.0);
        temperature_solver.Solve()
        print "after solving for the temperature"
        NistParameters.CalculateViscosity(fluid_model_part.Nodes) #applying the viscosity

##        #clearing the model part
        (temperature_model_part.Nodes).clear()
        (temperature_model_part.Elements).clear()
        (temperature_model_part.Conditions).clear()

##        for node in fluid_model_part.Nodes:
##            if(node.GetSolutionStepValue(VISCOSITY) > 1.0):
##                node.Fix(DISPLACEMENT_X)
##                node.Fix(DISPLACEMENT_Y)
##                node.Fix(DISPLACEMENT_Z)
##            else:
##                if(node.GetSolutionStepValue(IS_STRUCTURE) != 1):
##                    node.Free(DISPLACEMENT_X)
##                    node.Free(DISPLACEMENT_Y)
##                    node.Free(DISPLACEMENT_Z)
                    
         #finalizing fluid calculation step
        solver.Solve()
        
        print "after completing the solution"

        PrintVolume (time, scalarout, volume0)
        PrintVolumeCatchpan (time, catchpanout, volume0)

        if(time > next_output_time and time > 150.0):
    
            file_name = "firstNISTtest_lagrangian"
            file_name = file_name + str(time)
        
            gid_io.ChangeOutputName(file_name,GiDPostMode.GiD_PostBinary);
            gid_io.WriteMesh2D((combined_model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);

            gid_io.WriteNodalResults(DISPLACEMENT, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(NODAL_H, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FLUID, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_BOUNDARY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FREE_SURFACE, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(IS_STRUCTURE, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(VELOCITY, combined_model_part.Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(BULK_MODULUS, (combined_model_part).Nodes, time, 0);

            gid_io.WriteNodalResults(NODAL_AREA, (combined_model_part).Nodes, time, 0);

            gid_io.WriteNodalResults(VISCOSITY, (combined_model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(TEMPERATURE, (combined_model_part).Nodes, time, 0);

            gid_io.Flush()
            #gid_io.CloseResultFile();

            next_output_time = next_output_time  + output_Dt;

    print "completed step ",step
 


    

##
##
##
##
##Dt = 0.005
##nsteps = 10000
###output_Dt = 0.05
##output_Dt = output_dt
##min_dt = 0.00001
##max_dt = max_delta_time
##safety_factor = 0.5 #you should put a safety factor ;-)!!!
##
##next_output_time = output_Dt
##
###initializing the solver
##solver.Initialize()
##
##time = 0.0
##step = 0
##
##
##while (time < total_time):
##    step = step+1
##    
##    
##    
##    
##    print time
##    if(step <= 3):
##        new_Dt = 0.00000001;
##        time = time + new_Dt*safety_factor
##
##    #solving the fluid problem
##    if(step > 3):
##        new_Dt = solver.EstimateDeltaTime(max_dt,domain_size)
##        time = time + new_Dt*safety_factor
##
##        fluid_model_part.CloneTimeStep(time)
##        structure_model_part.CloneTimeStep(time)
##        combined_model_part.CloneTimeStep(time)
##
##        print "before the solution"
##
##        solver.Solve()
##        
##        print "after completing the solution"
##
##        if(time > next_output_time):
##    
##            file_name = project_name
##            file_name = file_name + str(time)
##        
##            gid_io.ChangeOutputName(file_name,GiDPostMode.GiD_PostBinary);
##            gid_io.WriteMesh((combined_model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
##
##            gid_io.WriteNodalResults(DISPLACEMENT, combined_model_part.Nodes, time, 0);
##            gid_io.WriteNodalResults(NODAL_H, combined_model_part.Nodes, time, 0);
##            gid_io.WriteNodalResults(IS_FLUID, combined_model_part.Nodes, time, 0);
##            gid_io.WriteNodalResults(IS_BOUNDARY, combined_model_part.Nodes, time, 0);
##            gid_io.WriteNodalResults(IS_FREE_SURFACE, combined_model_part.Nodes, time, 0);
##            gid_io.WriteNodalResults(IS_STRUCTURE, combined_model_part.Nodes, time, 0);
##            gid_io.WriteNodalResults(VELOCITY, combined_model_part.Nodes, time, 0);
##            gid_io.WriteNodalResults(PRESSURE, (combined_model_part).Nodes, time, 0);
##            gid_io.WriteNodalResults(BODY_FORCE, (combined_model_part).Nodes, time, 0);
##
##            gid_io.Flush()
##            #gid_io.CloseResultFile();
##
##            next_output_time = next_output_time  + output_Dt;
## 
##
##          
##        
##
##
