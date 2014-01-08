import time as timer
import os
import sys
import math
from numpy import *

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import datetime

import DEM_explicit_solver_var as DEM_parameters
import DEM_procedures
proc = DEM_procedures.Procedures(DEM_parameters)

import pressure_script as Press

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

# Defining a model part for the solid part

my_timer = Timer();
balls_model_part = ModelPart("SolidPart");

RigidFace_model_part   = ModelPart("RigidFace_Part");  
mixed_model_part       = ModelPart("Mixed_Part");

RigidFace_model_part.AddNodalSolutionStepVariable(VELOCITY)
RigidFace_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)

# Importing the strategy object

import continuum_sphere_strategy as SolverStrategy

SolverStrategy.AddVariables(balls_model_part, DEM_parameters)

# Reading the model_part: binary or ascii, multifile or single --> only binary and single for mpi.

if (DEM_parameters.OutputFileType == "Binary"):
    gid_mode = GiDPostMode.GiD_PostBinary
  
else:
    gid_mode = GiDPostMode.GiD_PostAscii
  
if (DEM_parameters.Multifile == "multiple_files"):
    multifile = MultiFileFlag.MultipleFiles
  
else:
    multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions   = WriteConditionsFlag.WriteConditions

gid_io = GidIO(DEM_parameters.problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(DEM_parameters.problem_name)
model_part_io_solid.ReadModelPart(balls_model_part)

rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
model_part_io_solid = ModelPartIO(rigidFace_mp_filename)
model_part_io_solid.ReadModelPart(RigidFace_model_part)

# Setting up the buffer size: SHOULD BE DONE AFTER READING!!!

balls_model_part.SetBufferSize(1)

# Adding dofs

SolverStrategy.AddDofs(balls_model_part)

# Constructing a creator/destructor object

creator_destructor = ParticleCreatorDestructor()

# Creating a solver object

solver = SolverStrategy.ExplicitStrategy(balls_model_part, RigidFace_model_part,creator_destructor, DEM_parameters) #here, solver variables initialize as default

# Creating necessary directories

main_path        = os.getcwd()
post_path        = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Post_Files'
list_path        = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Post_Lists'
neigh_list_path  = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Neigh_Lists'
data_and_results = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Results_and_Data'
graphs_path      = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Graphs'
MPI_results      = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_MPI_results'

for directory in [post_path, list_path, neigh_list_path, data_and_results, graphs_path, MPI_results]:

    if not os.path.isdir(directory):
        os.makedirs(str(directory))

os.chdir(list_path)

multifile        = open(DEM_parameters.problem_name + '_all' + '.post.lst', 'w')
multifile_5      = open(DEM_parameters.problem_name + '_5'   + '.post.lst', 'w')
multifile_10     = open(DEM_parameters.problem_name + '_10'  + '.post.lst', 'w')
multifile_50     = open(DEM_parameters.problem_name + '_50'  + '.post.lst', 'w')

multifile.write('Multiple\n')
multifile_5.write('Multiple\n')
multifile_10.write('Multiple\n')
multifile_50.write('Multiple\n')

first_print     = True
index_5         = 1
index_10        = 1
index_50        = 1
prev_time       = 0.0
control         = 0.0


os.chdir(main_path)

print 'Initializing Problem....'

solver.Initialize()

if ( (DEM_parameters.ContinuumOption =="ON") and (DEM_parameters.ContactMeshOption =="ON") ) :

  contact_model_part = solver.contact_model_part   

#-------------------------------------------------------------------------------------------------------------------------

#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALITZATIONS--------------------------------------------------------

Pressure = DEM_parameters.ConfinementPressure * 1e6 #Mpa

if (DEM_parameters.PredefinedSkinOption == "ON" ):

   proc.SetPredefinedSkin(balls_model_part)


if ( (DEM_parameters.ContinuumOption == "ON") and (DEM_parameters.ConcreteTestOption != "OFF") ):
  
    (sup_layer_fm, inf_layer_fm, sup_plate_fm, inf_plate_fm) = proc.ListDefinition(balls_model_part,solver)  # defines the lists where we measure forces

    strain=0.0; total_stress = 0.0; first_time_entry = 1
    # for the graph plotting    
    velocity_node_y = 0.0
    height = DEM_parameters.SpecimenHeight
    diameter = DEM_parameters.SpecimenWidth

    initial_time = datetime.datetime.now()

    os.chdir(graphs_path)

    chart = open("Provisional_CHART.grf", 'w')
    
    if(DEM_parameters.ConcreteTestOption == "BTS"):

        bts_export = open(DEM_parameters.problem_name +"_bts_"+str(datetime.datetime.now())+".grf",'w');      
        proc.BtsSkinDetermination(balls_model_part,solver,DEM_parameters)
        
    else:

      graph_export_top = open(DEM_parameters.problem_name +"_Provisional_TOP.grf", 'w')
      graph_export_bot = open(DEM_parameters.problem_name +"_Provisional_BOT.grf", 'w')
      graph_export_mean = open(DEM_parameters.problem_name +"_Provisional_MEAN.grf", 'w')

      
      #measuring height:
      pre_utilities = PreUtilities(balls_model_part)
    
      (subtotal_top,weight_top) = pre_utilities.MeasureTopHeigh(balls_model_part)
      (subtotal_bot,weight_bot) = pre_utilities.MeasureBotHeigh(balls_model_part)

      mean_top = subtotal_top/weight_top;
      mean_bot = subtotal_bot/weight_bot;
    
      ini_height = mean_top - mean_bot

      height = ini_height    
    
      print ('Initial Height of the Model: ' + str(ini_height)+'\n')
      
      
      if(DEM_parameters.PredefinedSkinOption == "ON" ):
        print "ERROR: in Concrete Test Option the Skin is automatically predefined. Switch the Predefined Skin Option OFF"

      (xtop_area,xbot_area,xlat_area,xtopcorner_area,xbotcorner_area) = proc.CylinderSkinDetermination(balls_model_part,solver,DEM_parameters) # defines the skin and areas
       
      if (DEM_parameters.PoissonMeasure =="ON"):
        
        graph_export_poisson = open(DEM_parameters.problem_name + "_poisson_"+str(datetime.datetime.now())+".csv",'w');

    os.chdir(main_path)
      
    
    if ( DEM_parameters.ConcreteTestOption == "TRIAXIAL" ) and (Pressure != 0.0) :

      #Correction Coefs
      alpha_top = 3.141592*diameter*diameter*0.25/(xtop_area + 0.70710678*xtopcorner_area)
      alpha_bot = 3.141592*diameter*diameter*0.25/(xbot_area + 0.70710678*xbotcorner_area)
      alpha_lat = 3.141592*diameter*height/(xlat_area + 0.70710678*xtopcorner_area + 0.70710678*xbotcorner_area) 

      print "Applying Pressure" , "\n"
 
      Press.ApplyPressure(Pressure, proc.XLAT, proc.XBOT, proc.XTOP, proc.XBOTCORNER, proc.XTOPCORNER,alpha_top,alpha_bot,alpha_lat)
      renew_pressure = 0
    
#
    


# Initialization of physics monitor and of the initial position of the center of mass

#physics_calculator = SphericElementGlobalPhysicsCalculator(balls_model_part)
#properties_list = []

print 'Initialitzation Complete' + '\n'


step                   = 0
time                   = 0.0
time_old_print         = 0.0
initial_pr_time        = timer.clock()
initial_real_time      = timer.time()

#-------------------------------------------------------------------------------------------------------------------------------------

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------

post_utility = PostUtilities()

os.chdir(post_path)

if (DEM_parameters.Multifile == "single_file"):

  post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
  if (DEM_parameters.ContactMeshOption == "ON"):
      post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)
  post_utility.AddModelPartToModelPart(mixed_model_part, RigidFace_model_part)
  gid_io.InitializeMesh(0.0) 
  gid_io.WriteMesh(RigidFace_model_part.GetMesh())
  gid_io.WriteSphereMesh(balls_model_part.GetMesh())
  if (DEM_parameters.ContactMeshOption == "ON"):
      gid_io.WriteMesh(contact_model_part.GetMesh())
  gid_io.FinalizeMesh()
  gid_io.InitializeResults(0.0, mixed_model_part.GetMesh())

  
##OEDOMETRIC

if(DEM_parameters.ConcreteTestOption == "OEDOMETRIC"):
  
  for node in proc.LAT:

    node.SetSolutionStepValue(VELOCITY_X,0.0);
    node.SetSolutionStepValue(VELOCITY_Z,0.0);
    node.Fix(VELOCITY_X);
    node.Fix(VELOCITY_Z);
    
##MODEL DATA 

if (DEM_parameters.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    (coordination_number) = proc.ModelData(balls_model_part, contact_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
    print ('Coordination Number: ' + str(coordination_number) + '\n' )
    os.chdir(main_path)

  
##WHEATHERFORD

w_densi = DEM_parameters.w_densi
w_dynfrc = DEM_parameters.w_dynfrc
w_young = DEM_parameters.w_young
w_poiss = DEM_parameters.w_poiss

os.chdir(graphs_path)

if(DEM_parameters.Dempack and (DEM_parameters.ConcreteTestOption != "OFF")):
  
  print("This chart is only valid for one material case" + '\n')

  chart.write(("***********PARAMETERS*****************")+'\n')
  chart.write( "                                    *" +'\n')
  chart.write( "    DENSI  = " + (str(w_densi))+" Kg/m3     "+'\n')
  chart.write( "    STAFRC = " + (str(DEM_parameters.InternalFriction))+"           "+'\n')
  chart.write( "    DYNFRC = " + (str(w_dynfrc))+"          " +'\n')
  chart.write( "    YOUNG  = " + (str(w_young/1e9))+" GPa"+"     " +'\n')
  chart.write( "    POISS  = " + (str(w_poiss))+"           " +'\n')
  chart.write( "    NTSR   = " + (str(DEM_parameters.SigmaMin))+" Mpa        " +'\n')
  chart.write( "    LCS1   = " + (str(DEM_parameters.C1))+" Mpa       " +'\n')
  chart.write( "    LCS2   = " + (str(DEM_parameters.C2))+" Mpa       " +'\n')
  chart.write( "    LCS3   = " + (str(DEM_parameters.C3))+" Mpa       " +'\n')
  chart.write( "    YRC1   = " + (str(DEM_parameters.N1))+"           " +'\n')
  chart.write( "    YRC2   = " + (str(DEM_parameters.N2))+"           " +'\n')
  chart.write( "    YRC3   = " + (str(DEM_parameters.N3))+"           " +'\n')
  chart.write( "    NG     = " + (str(7.0/6.0*2.0*(1.0+w_poiss)))+"           " +'\n')
  chart.write( "    FSS    = " + (str(DEM_parameters.TauZero))+" Mpa       " +'\n')
  chart.write( "    YEP    = " + (str(DEM_parameters.PlasticYoungModulus/1e9))+" GPa"+"     " +'\n')
  chart.write( "    YIELD  = " + (str(DEM_parameters.PlasticYieldStress))+" Mpa       " +'\n')
  chart.write( "    EDR    = " + (str(DEM_parameters.DamageDeformationFactor))+"           " +'\n')
  chart.write( "    GDAMP  = " + (str(DEM_parameters.DempackGlobalDamping))+"           " +'\n')
  chart.write( "    LDAMP  = " + (str(DEM_parameters.DempackDamping))+"           " +'\n')
  chart.write( "    ALPHA  = " + (str(DEM_parameters.AreaFactor))+"           " +'\n')
  chart.write( "                                    *" +'\n')
  chart.write( "**************************************" +'\n')

  chart.close()
  a_chart = open("Provisional_CHART.grf","r")
  
  for line in a_chart.readlines():
    print(line)
  a_chart.close()
  

#------------------------------------------------------------------------------------------
 
###########################################################################################
#                                                                                         #
#                                    MAIN LOOP                                            #
#                                                                                         #
###########################################################################################
os.chdir(main_path)

dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME

total_steps_expected = int(DEM_parameters.FinalTime / dt)

print ('Main loop starts at instant: ' + str(initial_pr_time) + '\n')

print ('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' + '\n' )

left_nodes = list()
right_nodes = list()

xleft_weight  = 0.0         
xright_weight  = 0.0

left_counter = 0.0
right_counter = 0.0

if(DEM_parameters.PoissonMeasure == "ON"):
        
      for node in balls_model_part.Nodes:
        
        if (node.GetSolutionStepValue(GROUP_ID)==4):
          
           left_nodes.append(node)
           xleft_weight = +(node.X0 - node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
           left_counter = +node.GetSolutionStepValue(RADIUS)
           
        elif(node.GetSolutionStepValue(GROUP_ID)==8):
          
           right_nodes.append(node)
           xright_weight = +(node.X + node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
           right_counter = +node.GetSolutionStepValue(RADIUS)
           
      width_ini = xright_weight/right_counter - xleft_weight/left_counter
        
step_to_fix_velocities = balls_model_part.ProcessInfo[STEP_FIX_VELOCITIES]
        
while (time < DEM_parameters.FinalTime):

    dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time = time + dt
    #balls_model_part.CloneTimeStep(time)
    balls_model_part.ProcessInfo[TIME] = time
    balls_model_part.ProcessInfo[DELTA_TIME] = dt
    balls_model_part.ProcessInfo[TIME_STEPS] = step

    #########################_SOLVE_#########################################4
    os.chdir(main_path)
    solver.Solve()
    #########################TIME CONTROL######################################4

    incremental_time = (timer.time() - initial_real_time) - prev_time

    if (incremental_time > DEM_parameters.ControlTime):
        percentage = 100 * (float(step) / total_steps_expected)

        print 'Real time calculation: ' + str(timer.time() - initial_real_time)
        print 'Percentage Completed: '  + str(percentage) + ' %'
        print "TIME STEP = "            + str(step) + '\n'
        
        if( DEM_parameters.ContinuumOption =="ON" and ( step >= step_to_fix_velocities ) and DEM_parameters.ConcreteTestOption != "OFF" and DEM_parameters.MonitoringOption == "ON"):        
            monitoring = PostUtilities().QuasiStaticAdimensionalNumber(balls_model_part,contact_model_part,balls_model_part.ProcessInfo)
            print "The quasi-static-adimensional-number is:  "            + str(monitoring) + '\n'
            print "The measured stiffness is:  "            + str(total_stress/strain/1e6) + "Mpa" + '\n'
            
        sys.stdout.flush()
        
        prev_time = (timer.time() - initial_real_time)
  
    if ((timer.time() - initial_real_time > 60) and first_print == True and step != 0):    
        first_print = False    
        estimated_sim_duration = 60.0 * (total_steps_expected / step) # seconds

        print('The calculation total estimated time is ' + str(estimated_sim_duration) + 'seconds' + '\n')
        print('in minutes:'        + str(estimated_sim_duration / 60.0) + 'min.' + '\n')
        print('in hours:'        + str(estimated_sim_duration / 3600.0) + 'hrs.' + '\n')
        print('in days:'        + str(estimated_sim_duration / 86400.0) + 'days' + '\n') 
        sys.stdout.flush()
        
        if (estimated_sim_duration / 86400 > 2.0):
            print('WARNING!!!:       VERY LASTING CALCULATION' + '\n')

    #########################CONCRETE_TEST_STUFF#########################################4

    os.chdir(data_and_results)
                                                                                                                                                                                               
    if( (DEM_parameters.ConcreteTestOption == "TRIAXIAL" ) and (Pressure != 0.0) ):
       #and (step < 0.01*DEM_parameters.TotalTimePercentAsForceAplTime*total_steps_expected) )
        
        if( renew_pressure == 10):
          
          Press.ApplyPressure(Pressure, proc.XLAT, proc.XBOT, proc.XTOP, proc.XBOTCORNER, proc.XTOPCORNER,alpha_top,alpha_bot,alpha_lat)
                 
          renew_pressure = 0
    
        renew_pressure += 1
    
    total_force_top = 0.0
    total_force_bot = 0.0
    total_force_bts = 0.0
    
    if( DEM_parameters.ContinuumOption =="ON"):
      
        if( DEM_parameters.ConcreteTestOption =="BTS"):

          for node in sup_layer_fm:

            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]

            total_force_bts += force_node_y
        
        elif ( ( step >= step_to_fix_velocities ) and DEM_parameters.ConcreteTestOption != "OFF"):

          if(first_time_entry):
            #measuring height:

            pre_utilities = PreUtilities(balls_model_part)

            (subtotal_top,weight_top) = pre_utilities.MeasureTopHeigh(balls_model_part)
            (subtotal_bot,weight_bot) = pre_utilities.MeasureBotHeigh(balls_model_part)

            mean_top = subtotal_top/weight_top;
            mean_bot = subtotal_bot/weight_bot;

            ini_height2 = mean_top - mean_bot

            print 'Current Height after confinement: ' + str(ini_height2) + '\n'
            print 'Axial strain due to the confinement: ' + str( 100*(ini_height2-ini_height)/ini_height ) + ' %' +'\n'
            height = ini_height2
          
            for node in sup_layer_fm:
              velocity_node_y = node.GetSolutionStepValue(VELOCITY_Y) #Applied velocity during the uniaxial compression test
              break

            print 'velocity for the graph: ' + str(velocity_node_y) + '\n'
                  
            first_time_entry = 0

          strain += -1.0*velocity_node_y*dt/height

          for node in sup_layer_fm:

            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]

            
            total_force_top += force_node_y

          total_stress_top = total_force_top/(DEM_parameters.MeasuringSurface*1000000)
          
          for node in inf_layer_fm:

            force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]

            total_force_bot += force_node_y

          total_stress_bot = total_force_bot/(DEM_parameters.MeasuringSurface*1000000)
          
          if(DEM_parameters.PoissonMeasure == "ON"):
                      
            xleft_weight  = 0.0         
            xright_weight  = 0.0

            left_counter = 0.0
            right_counter = 0.0

            for node in left_nodes:
              
              xleft_weight = +(node.X - node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
              left_counter = +node.GetSolutionStepValue(RADIUS)
              
            for node in right_nodes:
              
              xright_weight = +(node.X + node.GetSolutionStepValue(RADIUS))*node.GetSolutionStepValue(RADIUS)
              right_counter = +node.GetSolutionStepValue(RADIUS)
            
            width_now = xright_weight/right_counter - xleft_weight/left_counter

            measured_poisson =  ((width_now-width_ini)/width_ini)/strain

            #print( (width_now/0.05)/strain )

    os.chdir(list_path)    
    multifile.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')   
    os.chdir(main_path)

  ##########################___GiD IO____#########################################4

    time_to_print = time - time_old_print

    if (time_to_print >= DEM_parameters.OutputTimeStep):
        os.chdir(data_and_results)

        #properties_list = proc.MonitorPhysicalProperties(balls_model_part, physics_calculator, properties_list)

        if (index_5 == 5):
            multifile_5.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_5 = 0

        if (index_10 == 10):
            multifile_10.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_10 = 0

        if (index_50 == 50):
            multifile_50.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_50 = 0

        index_5  += 1
        index_10 += 1
        index_50 += 1

        if (DEM_parameters.Multifile == "multiple_files"):
            gid_io.FinalizeResults()

        os.chdir(graphs_path)

        if (DEM_parameters.PrintNeighbourLists == "ON"): # Printing neighbours id's
            os.chdir(neigh_list_path)
            neighbours_list = open('neigh_list_' + str(time), 'w')

            for elem in balls_model_part.Elements:
                ID = (elem.Id)
                Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)

            for i in range(len(Neigh_ID)):
                neighbours_list.write(str(ID) + ' ' + str(Neigh_ID[i]) + '\n')

            neighbours_list.close()

        os.chdir(post_path)

        if (DEM_parameters.Multifile == "multiple_files"):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()

            post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
            if (DEM_parameters.ContactMeshOption == "ON"):
                post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)
            post_utility.AddModelPartToModelPart(mixed_model_part, RigidFace_model_part)

            gid_io.InitializeMesh(time) 
            gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            if (DEM_parameters.ContactMeshOption == "ON"):
                gid_io.WriteMesh(contact_model_part.GetMesh())
            gid_io.WriteMesh(RigidFace_model_part.GetMesh())
            gid_io.FinalizeMesh()
            
            gid_io.InitializeResults(time, mixed_model_part.GetMesh())
                                
        proc.PrintingGlobalVariables(gid_io, mixed_model_part, time)
        proc.PrintingBallsVariables(gid_io, balls_model_part, time)
        if (DEM_parameters.ContactMeshOption == "ON"):
            proc.PrintingContactElementsVariables(gid_io, contact_model_part, time)
        
        os.chdir(main_path)     
              
        if (DEM_parameters.Multifile == "multiple_files"):
            gid_io.FinalizeResults() 

        time_old_print = time
    
    if(DEM_parameters.ConcreteTestOption == "BTS"):
       bts_export.write(str(step)+"  "+str(total_force_bts)+'\n')
      
    if ( ( (DEM_parameters.ConcreteTestOption =="TRIAXIAL") or (DEM_parameters.ConcreteTestOption == "UCS")) and (step >= step_to_fix_velocities )):
      graph_export_top.write(str(strain)+"  "+str(total_stress_top)+'\n')
      graph_export_bot.write(str(strain)+"  "+str(total_stress_bot)+'\n')
      total_stress_mean = 0.5*(total_stress_bot + total_stress_top)
      graph_export_mean.write(str(strain)+"  "+str(total_stress_mean)+'\n')
      
      if (DEM_parameters.PoissonMeasure =="ON"):
        graph_export_poisson.write(str(strain)+"  "+str(measured_poisson)+'\n')
      
         
    step += 1
#-------------------------------------------------------------------------------------------------------------------------------------


#-----------------------FINALITZATION OPERATIONS-------------------------------------------------------------------------------------- 
#proc.PlotPhysicalProperties(properties_list, graphs_path)

if (DEM_parameters.Multifile == "single_file"):
    gid_io.FinalizeResults()


os.chdir(graphs_path)

#for filename in os.listdir("."):
#  if filename.startswith("Provisional_TOP"):
#    os.rename(filename, DEM_parameters.problem_name + "_graph_" + str(initial_time) + "_TOP.grf")
#  if filename.startswith("Provisional_BOT"):
#    os.rename(filename, DEM_parameters.problem_name + "_graph_" + str(initial_time) + "_BOT.grf")
#  if filename.startswith("Provisional_MEAN"):
#    os.rename(filename, DEM_parameters.problem_name + "_graph_" + str(initial_time) + "_MEAN.grf")
#  if filename.startswith("Provisional_CHART"):
#    os.rename(filename, DEM_parameters.problem_name + "_CHART_" + str(initial_time) +".grf")

for filename in os.listdir("."):
  if filename.startswith(DEM_parameters.problem_name +"_Provisional_TOP"):
    os.rename(filename, DEM_parameters.problem_name + "_graph_" + "_TOP.grf")
  if filename.startswith(DEM_parameters.problem_name +"_Provisional_BOT"):                    
    os.rename(filename, DEM_parameters.problem_name + "_graph_" + "_BOT.grf")
  if filename.startswith(DEM_parameters.problem_name +"_Provisional_MEAN"):                   
    os.rename(filename, DEM_parameters.problem_name + "_graph_" + "_MEAN.grf")
  if filename.startswith(DEM_parameters.problem_name +"_Provisional_CHART"):                  
    os.rename(filename, DEM_parameters.problem_name + "_CHART_" +".grf")


if (DEM_parameters.ConcreteTestOption!= "OFF"):
  graph_export_top.close()
  graph_export_bot.close()
  graph_export_mean.close()

  
if(DEM_parameters.ConcreteTestOption == "BTS"):
  bts_export.close()
  
multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()
os.chdir(main_path)

elapsed_pr_time     = timer.clock() - initial_pr_time
elapsed_real_time   = timer.time() - initial_real_time

print 'Calculation ends at instant: '                 + str(timer.time())
print 'Calculation ends at processing time instant: ' + str(timer.clock())
print 'Elapsed processing time: '                     + str(elapsed_pr_time)
print 'Elapsed real time: '                           + str(elapsed_real_time)

print (my_timer)

print "ANALYSIS COMPLETED"
