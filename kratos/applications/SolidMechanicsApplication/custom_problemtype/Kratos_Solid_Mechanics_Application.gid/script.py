#Import the general variables read from the GiD
import problem_settings as general_variables

# Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/jmaria/kratos')


##################################################################
### ***************GENERAL MAIN OF THE ANALISYS****************###
##################################################################

##time control starts
from time import *
print ctime()
# measure process time
t0p = clock()
# measure wall time
#t0w = time()

######################--CONFIGURATIONS START--####################
#setting the domain size for the problem to be solved
domain_size = general_variables.domain_size

#including kratos path
from KratosMultiphysics import *

#including Applications path
from KratosMultiphysics.SolidMechanicsApplication import *
if(general_variables.LinearSolver != "SkylineLUFactorization" and general_variables.LinearSolver != "ParallelMKLPardisoSolver" ): 
  from KratosMultiphysics.ExternalSolversApplication import *

if(general_variables.LinearSolver == "ParallelMKLPardisoSolver"):
  from KratosMultiphysics.MKLSolversApplication import *
  

######################--SET NUMBER OF THREADS --#########
#if the data units are expressed in the INTERNATIONAL system
def SetParallelSize(num_threads):
  parallel=OpenMPUtils()
  parallel.SetNumThreads(num_threads); 
######################--SET NUMBER OF THREADS --############

#defining the number of threads:
num_threads = general_variables.NumberofThreads
SetParallelSize(num_threads)


#defining the name of the problem:
problem_name = general_variables.problem_name
problem_path = general_variables.problem_path

#defining a model part
model_part = ModelPart("SolidDomain");  

if(general_variables.LoadRestart == "True"):
  restart_path= problem_path + "/" + problem_name + "_" + str(general_variables.Restart_Step)
  #serializer = Serializer(restart_path,SERIALIZER_TRACE_ERROR);
  serializer = Serializer(restart_path,SERIALIZER_NO_TRACE);
  serializer.Load("ModelPart",model_part);
else:
  #add specific variables for the problem (conditions) 
  model_part.AddNodalSolutionStepVariable(IMPOSED_DISPLACEMENT);
  model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
  model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
  model_part.AddNodalSolutionStepVariable(FORCE);
  model_part.AddNodalSolutionStepVariable(FACE_LOAD);
  model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

  
#creating auxiliar list for printing list files:
initial_list  = True
linitial_list = []
listprint = []  
for lfile in general_variables.file_list:
  listprint.append(1);
  linitial_list.append(True);


if(general_variables.LoadRestart == "False"):
  print "Remove ALL previous FILES"
  #remove previous results:
  filelist1 = [ f for f in os.listdir(problem_path) if f.endswith(".bin") ]
  for f in filelist1:
    os.remove(f)

  #remove previous list files:
  filelist2 = [ f for f in os.listdir(problem_path) if f.endswith(".lst") ]
  for f in filelist2:
    os.remove(f)

  #remove previous restart files:
  filelist3 = [ f for f in os.listdir(problem_path) if f.endswith(".rest") ]
  for f in filelist3:
    os.remove(f)

  #remove previous graph files:
  filelist4 = [ f for f in os.listdir(problem_path) if f.endswith(".png") ]
  for f in filelist4:
    os.remove(f)
else:
  print "Remove RESTART posterior FILES"
  #remove previous results after restart:
  total_files = 0;
  filelist1 = []
  filelist2 = []
  for f in os.listdir(problem_path):
    if(f.endswith(".post.bin")):
      total_files = total_files + 1
     
  total_files = total_files - general_variables.Restart_Step + 250  #arbitrary number to ensure to remove all rest files

  for rfile in range(0,total_files):
    for f in os.listdir(problem_path):
      if(f.endswith("_"+str(general_variables.Restart_Step+rfile+1)+".post.bin")):
        filelist1.append(f);
      if(f.endswith("_"+str(general_variables.Restart_Step+rfile+1)+".png")):
        filelist2.append(f);
 
      
  for f in filelist1:
    os.remove(f)

  for f in filelist2:
    os.remove(f)


  #remove previous restart files:
  counter = 1;
  total_files = 0;
  filelist2 = []
  for f in os.listdir(problem_path):
    if(f.endswith(".post.bin")):
      total_files = total_files + 1
     
  total_files = total_files + 100  #arbitrary number to ensure to remove all rest files
      
  for rfile in range(0,total_files):
    for f in os.listdir(problem_path):
      if(f.endswith("_"+str(general_variables.Restart_Step+rfile+1)+".rest")):
        filelist2.append(f);
 
  for f in filelist2:
    os.remove(f)

  #remove previous list files and rebuild it:
  filelist2 = [ f for f in os.listdir(problem_path) if f.endswith(".post.lst") ]
  for f in filelist2:
    os.remove(f)
    
  #print list files:
  if(general_variables.PrintLists == "True"):
    total_files = 0;
    for f in os.listdir(problem_path):
      if(f.endswith(".post.bin")):
        total_files = total_files + 1
     
    for rfile in range(0,total_files):
      for f in os.listdir(problem_path):
        if(f.endswith("_"+str(rfile)+".post.bin")):
          problempath= problem_path + "/" + problem_name + "_1.post.lst"
          if(os.path.exists(problempath) == False):
            listfile = open(problempath,"a")
            problemname = "Multiple\n"
            listfile.write(problemname)
            problemname = problem_name + "_" + str(rfile) + ".post.bin\n" 
            listfile.write(problemname)
            listfile.close()
            file_set = True
          else:
            listfile = open(problempath,"a")
            problemname = problem_name + "_" + str(rfile) + ".post.bin\n" 
            listfile.write(problemname)
            listfile.close()
            file_set = True
            
          num_list_files = len(general_variables.file_list) 
          for lfile in range(0,num_list_files):
            if( general_variables.file_list[lfile] == listprint[lfile] ):
              problempath= problem_path + "/" + problem_name + "_"+ str(general_variables.file_list[lfile]) + ".post.lst"
              if(os.path.exists(problempath) == False):
                listfile = open(problempath,"a")
                problemname = "Multiple\n" 
                listfile.write(problemname)
                linitial_list[lfile] = False
                problemname = problem_name + "_" + str(rfile) + ".post.bin\n" 
                listfile.write(problemname)
                listfile.close()
                listprint[lfile] = 1
              else:
                listfile = open(problempath,"a")
                problemname = problem_name + "_" + str(rfile) + ".post.bin\n" 
                listfile.write(problemname)
                listfile.close()
                listprint[lfile] = 1
            else:
              listprint[lfile] = listprint[lfile]+1
        
  

######################--CONFIGURATIONS END--####################


######################--SPECIFIC SOLVER START--#################
#x = raw_input("stopped to allow debug: set breakpoints and press enter to continue");

#including phyton solver type
if(general_variables.SolverType == "StaticSolver"): 
  print " Static Solver "
  import solid_mechanics_main_solver as general_solver
  dynamicstype = 0; 
elif(general_variables.SolverType == "DynamicSolver"):
  print " Dynamic Solver "
  import solid_mechanics_main_solver as general_solver
  dynamicstype = 1;

#including solver type variables
if(general_variables.LoadRestart == "False"):
  general_solver.AddVariables(model_part)   

######################--SPECIFIC SOLVER END--####################


######################--GID OUTPUT OPTIONS START--###############
#setting the output on GID:
#if(general_variables.gid_output_mode == "Binary"):
gid_mode = GiDPostMode.GiD_PostBinary
#else:
#  gid_mode = GiDPostMode.GiD_PostAscii

gid_multiple_files  =  MultiFileFlag.MultipleFiles

gid_configuration   =  WriteDeformedMeshFlag.WriteDeformed
if(general_variables.WriteMesh == "Undeformed"):
  gid_configuration   =  WriteDeformedMeshFlag.WriteUndeformed

gid_conditions      =  WriteConditionsFlag.WriteElementsOnly
if(general_variables.WriteConditions == "True"):
  gid_conditions      =  WriteConditionsFlag.WriteConditions


gid_io = GidIO(problem_name, gid_mode, gid_multiple_files, gid_configuration, gid_conditions)
######################--GID OUTPUT OPTIONS END--#################


#------------------------#--FUNCTIONS START--#------------------#
#---------------------------------------------------------------#

######################--PRINTING TIME START --###################
def PrintSpentTime():
  print ctime()
  # measure process time
  tfp = clock()
  print "Analysis Spent Time  [Process Time = ",tfp - t0p,"] "
######################--PRINTING TIME END --###################


######################--CONSTITUTIVE LAW DEFINITION START--######
#set the constitutive law
def SetConstitutiveLaw():
  if(domain_size == 2):
    for prop in model_part.Properties:
      ConstitutiveLawName=prop.GetValue(CONSTITUTIVE_LAW_NAME)
      if(ConstitutiveLawName == "LinearElasticPlaneStrain"):
        prop.SetValue(CONSTITUTIVE_LAW, LinearElasticPlaneStrain2DLaw() )
        print "Linear Elastic Plane Strain 2D model selected"
      elif(ConstitutiveLawName == "LinearElasticPlaneStress"):
        prop.SetValue(CONSTITUTIVE_LAW, LinearElasticPlaneStress2DLaw() )
        print "Linear Elastic Plane Stress 2D model selected"
      elif(ConstitutiveLawName == "HyperElasticPlaneStrain"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticPlaneStrain2DLaw() )
        print "Hyperelastic Plane Strain 2D model selected"
      elif(ConstitutiveLawName == "HyperElasticPlaneStrainUP"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticUPPlaneStrain2DLaw() )
        print "Hyperelastic UP Plane Strain 2D model selected"
      elif(ConstitutiveLawName == "LinearElasticAxisym"):
        prop.SetValue(CONSTITUTIVE_LAW, LinearElasticAxisym2DLaw() )
        print "Linear Elastic Axisym 2D model selected"
      elif(ConstitutiveLawName == "HyperElasticAxisym"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticAxisym2DLaw() )
        print "Hyper Elastic Axisym 2D model selected"
      elif(ConstitutiveLawName == "HyperElasticAxisymUP"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticUPAxisym2DLaw() )
        print "Hyper Elastic UP Axisym 2D model selected"
      elif(ConstitutiveLawName == "HyperElasticPlasticJ2PlaneStrain"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticPlasticJ2PlaneStrain2DLaw() )
        print "Hyper Elastic Plastic J2 Plane Strain 2D model selected"
      else:
        print "ERROR: CONSTITUTIVE_LAW 2D not defined properly"
      
  else:
    for prop in model_part.Properties:
      ConstitutiveLawName=prop.GetValue(CONSTITUTIVE_LAW_NAME);
      if(ConstitutiveLawName == "LinearElastic"):
        prop.SetValue(CONSTITUTIVE_LAW, LinearElastic3DLaw() )
        print "Linear Elastic 3D model selected"
      elif(ConstitutiveLawName == "HyperElastic"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElastic3DLaw() )
        print "Hyperelastic 3D model selected"
      elif(ConstitutiveLawName == "HyperElasticUP"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticUP3DLaw() )
        print "Hyperelastic UP 3D model selected"
      elif(ConstitutiveLawName == "HyperElasticPlasticJ2"):
        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticPlasticJ23DLaw() )
        print "Hyper Elastic Plastic J2 3D model selected"
      else:
        print "ERROR: CONSTITUTIVE_LAW 3D not defined properly"



#set constitutive law (must be put in another python file)
##if(general_variables == "ThermoPlastic"):
##  model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, HyperElasticPlastic2D ())
######################--CONSTITUTIVE LAW DEFINITION END--########

######################--PROBLEM TYPE DEFINITION START--##########
#define the linear solver
def SetProblemType():
  if(general_variables.ProblemType == "Mechanical"):
    time_step_solver.ComputeMechanicsFlag = True
    time_step_solver.ComputeThermalFlag = False

######################--PROBLEM TYPE DEFINITION END--############

######################--LINEAR SOLVER DEFINITION START--#########
#define the linear solver
def SetLinearSolver():

  class linear_solver_config:
    solver_type = general_variables.LinearSolver
    scaling = False
    tolerance = general_variables.Linear_Solver_Tolerance #1e-7
    max_iteration = general_variables.Linear_Solver_Max_Iteration  #300     
    verbosity = 0
    is_symmetric = False  
    #pastix iterative
    gmres_krylov_space_dimension = 100
    ilu_level_of_fill = 3 #5
    #GMRES CG
    preconditioner_type = "None"
    #Deflated conjugate gradient
    assume_constant_structure = True
    max_reduced_size = 1000
    #AMG
    smoother_type = "ILU0" #"DAMPED_JACOBI"
    krylov_type = "GMRES"

  import linear_solver_factory

  time_step_solver = linear_solver_factory.ConstructSolver(linear_solver_config)
    
#Inside of the solver:
#builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)
#classical call
#self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,builder_and_solver,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   
#call passing builder_and_solver
#self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,builder_and_solver,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   


######################--LINEAR SOLVER DEFINITION END--##########


##########--CONVERGENCE CRITERION DEFINITION START--############
#DisplacementCriteria or ResidualCritera are from Kratos base
#DisplacementConvergenceCriteria or ResidualConvergenceCritera are from SolidMechanicsApplication
def SetConvergenceCriteria():
  #mechanical convergence criteria
  CT = general_variables.Convergence_Tolerance;
  AT = general_variables.Absolute_Tolerance; 
  if(general_variables.Convergence_Criteria == "Displacement_Criteria"):
    time_step_solver.mechanical_convergence_criteria  =  DisplacementConvergenceCriteria(CT,AT)
  elif(general_variables.Convergence_Criteria == "Residual_Criteria"): 
    time_step_solver.mechanical_convergence_criteria  =  ResidualConvergenceCriteria(CT,AT)
  elif(general_variables.Convergence_Criteria == "And_Criteria"): 
    Displacement   =   DisplacementConvergenceCriteria(CT,AT)
    Residual       =   ResidualConvergenceCriteria(CT,AT)
    time_step_solver.mechanical_convergence_criteria  = AndCriteria(Residual, Displacement)
  elif(general_variables.Convergence_Criteria == "Or_Criteria"): 
    Displacement   =   DisplacementConvergenceCriteria(CT,AT)
    Residual       =   ResidualConvergenceCriteria(CT,AT)
    time_step_solver.mechanical_convergence_criteria  = OrCriteria(Residual, Displacement)
  elif(general_variables.Convergence_Criteria == "Mixed_Criteria"):
    Displacement   =   MixedElementConvergenceCriteria(CT,AT)
    Residual       =   ResidualConvergenceCriteria(CT,AT)
    time_step_solver.mechanical_convergence_criteria  = AndCriteria(Residual, Displacement)
###########--CONVERGENCE CRITERION DEFINITION END--##############

###############--GID PRINT RESULTS DEFINITION START--############
def PrintResults(model_part,time,initial_time,step):
  #print "WRITING MODEL "
  print "WRITING RESULTS: [STEP: ",step,"] [TIME: ",time-initial_time,"]"
  gid_io.InitializeResults(step,(model_part).GetMesh())
  gid_io.InitializeMesh(step);
  if(general_variables.WriteParticles == "True"):
    gid_io.WriteNodeMesh((model_part).GetMesh());
  gid_io.WriteMesh((model_part).GetMesh());
  gid_io.FinalizeMesh();
  gid_time=time-initial_time
  if(general_variables.ProblemType == "Mechanical"):
    gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,gid_time,0)
    gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,gid_time,0)
    gid_io.WriteNodalResults(REACTION,model_part.Nodes,gid_time,0)
    gid_io.WriteNodalResults(FORCE_EXTERNAL,model_part.Nodes,gid_time,0)

    gid_io.PrintOnGaussPoints(CAUCHY_STRESS_TENSOR,model_part,gid_time)
    gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR,model_part,gid_time)
    gid_io.PrintOnGaussPoints(VON_MISES_STRESS,model_part,gid_time)
    gid_io.PrintOnGaussPoints(PLASTIC_STRAIN,model_part,gid_time)
    gid_io.PrintOnGaussPoints(DELTA_PLASTIC_STRAIN,model_part,gid_time)

    if(general_variables.SolverType == "DynamicSolver"):
      gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,gid_time,0)
      gid_io.WriteNodalResults(ACCELERATION,model_part.Nodes,gid_time,0)
  
  gid_io.Flush()
  gid_io.FinalizeResults()
  print " -Run_GID_for_viewing_results_of_the_analysis-"
  
  #print list files:
  if(general_variables.PrintLists == "True"):
    problempath= problem_path + "/" + problem_name + "_1.post.lst"
    if(os.path.exists(problempath) == False):
      listfile = open(problempath,"a")
    #if(step == 0 and general_variables.LoadRestart == "False"):
      problemname = "Multiple\n" 
      listfile.write(problemname)
    #endif
      problemname = problem_name + "_" + str(step) + ".post.bin\n" 
      listfile.write(problemname)
      listfile.close()
    else:
      listfile = open(problempath,"a")
      problemname = problem_name + "_" + str(step) + ".post.bin\n" 
      listfile.write(problemname)
      listfile.close()

    num_list_files = len(general_variables.file_list) 
    for lfile in range(0,num_list_files):
      if( general_variables.file_list[lfile] == listprint[lfile] ):
        problempath= problem_path + "/" + problem_name + "_"+ str(general_variables.file_list[lfile]) + ".post.lst"
        if(os.path.exists(problempath) == False):
          listfile = open(problempath,"a")
        #if(linitial_list[lfile] == True and general_variables.LoadRestart == "False"):
          problemname = "Multiple\n" 
          listfile.write(problemname)
          linitial_list[lfile] = False
        #endif
          problemname = problem_name + "_" + str(step) + ".post.bin\n" 
          listfile.write(problemname)
          listfile.close()
          listprint[lfile] = 1
        else:
          listfile = open(problempath,"a")
          problemname = problem_name + "_" + str(step) + ".post.bin\n" 
          listfile.write(problemname)
          listfile.close()
          listprint[lfile] = 1
      else:
        listprint[lfile] = listprint[lfile]+1



###############--GID PRINT RESULTS DEFINITION END--##############

#------------------------#--FUNCTIONS END--#--------------------#
#---------------------------------------------------------------#


######################--READ AND SET MODEL START--#######################
#define problem type
problemtype=general_variables.ProblemType
buffer_size=3;

if(general_variables.LoadRestart == "False"):
  #reading the model
  model_part_io = ModelPartIO(problem_name)
  model_part_io.ReadModelPart(model_part)

  #Note: the buffer size should be set up here
  # - once the mesh is read for the first time -

  #setting the buffer size
  model_part.SetBufferSize(buffer_size)

  #importing the solver files
  general_solver.AddDofs(model_part,problemtype)

  #writing the initial mesh
  #mesh_name = 0.0
  #gid_io.InitializeMesh( mesh_name );  
  #gid_io.WriteMesh((model_part).GetMesh());
  #gid_io.FinalizeMesh()

  #set the constitutive law
  SetConstitutiveLaw();


######################--READ AND SET MODEL END--########################



#--- PRINT CONTROL ---#
print model_part
print model_part.Properties


#########################--START SOLUTION--######################
#################################################################

#set time integration solver
time_step_solver = general_solver.SolidSolution(model_part,domain_size,problemtype,dynamicstype) 

#redefining the linear solver
SetProblemType();

#redefining the linear solver
SetLinearSolver();

#setting the convergence criteria
SetConvergenceCriteria();

#########################--INITIALIZE--###########################
##################################################################

## Assign variables to conditions
##for cond in model_part.Conditions:
##  cond.SetValue( PRESSURE, 1000 )

## Solver initialize
if(general_variables.LoadRestart == "True"):
  time_step_solver.Initialize(True)
else:
  time_step_solver.Initialize(False)

## Set variables
#set max iters
time_step_solver.SetMaxIters(int(general_variables.Max_Iter))
#set echo level
time_step_solver.SetEchoLevel(int(general_variables.echo_level));


#######################--TIME INTEGRATION--#######################
##################################################################
#define initial buffer needed steps
ini_steps = buffer_size-1  #to ensure the nodal variables initialization in buffer
#define loop range of steps
isteps    = 0
nsteps    = int(general_variables.nsteps) + ini_steps + 1 

#define time and time_step
time_step    =  general_variables.time_step
initial_time =  ini_steps*time_step;

if(general_variables.LoadRestart == "False"):
  model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;
  #incremental displacement
  if(general_variables.Incremental_Displacement == "True"):
    for node in model_part.Nodes:
      ImposedDisp  = node.GetSolutionStepValue(IMPOSED_DISPLACEMENT);
      Displacement = node.GetSolutionStepValue(DISPLACEMENT);
      Velocity     = node.GetSolutionStepValue(VELOCITY);
      #For displacement imposition
      if(node.IsFixed(DISPLACEMENT_X)==1):
        ImposedDisp[0]  = Displacement[0];
        Displacement[0] = 0;
      if(node.IsFixed(DISPLACEMENT_Y)==1):
        ImposedDisp[1]  = Displacement[1];
        Displacement[1] = 0;  
      if(node.IsFixed(DISPLACEMENT_Z)==1):
        ImposedDisp[2]  = Displacement[2];
        Displacement[2] = 0;  
      #For velocity imposition instead of displacement        
      if(node.IsFixed(VELOCITY_X)==1):
        ImposedDisp[0]  = Velocity[0]*time_step;
        Velocity[0] = 0;
      if(node.IsFixed(VELOCITY_Y)==1):
        ImposedDisp[1]  = Velocity[1]*time_step;
        Velocity[1] = 0;  
      if(node.IsFixed(VELOCITY_Z)==1):
        ImposedDisp[2]  = Velocity[2]*time_step;
        Velocity[2] = 0;   
      #print " ImposedDisp  =", ImposedDisp
      #print " Displacement =", Displacement
      #print " Velocity     =", Velocity
      node.SetSolutionStepValue(IMPOSED_DISPLACEMENT,ImposedDisp);
      #set to zero
      node.SetSolutionStepValue(DISPLACEMENT,Displacement);
      node.SetSolutionStepValue(VELOCITY,Velocity);

else:
  time_step =  model_part.ProcessInfo[DELTA_TIME];
  isteps    =  model_part.ProcessInfo[TIME_STEPS]+1;
  print "isteps: ",isteps


#writing and remeshing parameters
write_id    = 0;
count_steps = 0;
write_step  = 0;


#restart
step_to_restart  = 1
restart_interval = 10
if(general_variables.SaveRestart == "True"):
  restart_interval = general_variables.Restart_Interval;
if(general_variables.LoadRestart == "True"):
  step_to_restart = 0;
  count_steps = isteps - ini_steps;
  #recover meshing steps..etc
  write_id    = model_part.ProcessInfo[WRITE_ID];
  write_step  = count_steps;
  


for step in range(isteps,nsteps):
  time = time_step*step
  model_part.CloneTimeStep(time)
  model_part.ProcessInfo[DELTA_TIME] = time_step
  model_part.ProcessInfo[TIME_STEPS] = step

  if( (time-initial_time) > 0 ):
    print "STEP = ", step - ini_steps - 1
    print "TIME = ", time - initial_time
    model_part.ProcessInfo[TIME] = time - initial_time
  
  step_printed = False;

  #solving the solid problem
  if(step > ini_steps ):
    
    #solve time step non-linear system
    time_step_solver.Solve()

    print " STEP result stored "
      
    #print the results at the end of the step
    if(general_variables.WriteResults == "PreMeshing" and (count_steps == write_step or step == nsteps-1) ):
      #PrintResults(model_part, time, initial_time, int((step-ini_steps)));
      PrintResults(model_part, time, initial_time, int(write_id))
      if(write_step == 0):
        write_step = general_variables.WriteFrequency
      else:
        write_step = write_step+general_variables.WriteFrequency
      write_id = write_id+1
      step_printed = True;
      model_part.ProcessInfo[WRITE_ID] = write_id;
      PrintSpentTime();
    
   
    #update previous time step
    model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;
    
    #incremental load
    if(general_variables.Incremental_Load == "True"):
      for node in model_part.Nodes:
        force = node.GetSolutionStepValue(FACE_LOAD);
        force = force/(time_step*(step-ini_steps))
        force = force*time_step*(step-ini_steps+1)
        node.SetSolutionStepValue(FACE_LOAD,force);


    #print the results at the end of the remesh
    if(general_variables.WriteResults == "PostMeshing" and (count_steps == write_step or step == nsteps-1) ):
       #PrintResults(model_part, time, initial_time, int((step-ini_steps)));
      PrintResults(model_part, time, initial_time, int(write_id))
      if(write_step == 0):
        write_step = general_variables.WriteFrequency
      else:
        write_step = write_step+general_variables.WriteFrequency
      write_id   = write_id+1
      step_printed = True;
      model_part.ProcessInfo[WRITE_ID] = write_id;
      PrintSpentTime();
      

    #print restart file
    if( step_printed == True ):
      if( general_variables.SaveRestart == "True" and step_to_restart == restart_interval ): 
        restart_path= problem_path + "/" + problem_name + "_" + str(write_id-1)
        #serializer = Serializer(restart_path,SERIALIZER_TRACE_ALL);
        #serializer = Serializer(restart_path,SERIALIZER_TRACE_ERROR);
        serializer = Serializer(restart_path,SERIALIZER_NO_TRACE);        
        serializer.Save("ModelPart",model_part);
        serializer = 0;
        print " RESTART PRINTED :", write_id-1, " (Step: ", step, ") ";
        step_to_restart = 1;
      else:
        step_to_restart = step_to_restart + 1;

    count_steps = count_steps + 1
    
##########################--FINALIZE--############################
##################################################################

print "Analysis Finalized "
gid_io.FinalizeResults()

############################--END--###############################
##################################################################

# measure process time
tfp = clock()
# measure wall time
#tfw = time()

print ctime()
#print "Analysis Completed  [Process Time = ", tfp - t0p, "seconds, Wall Time = ", tfw - t0w, " ]"
print "Analysis Completed  [Process Time = ",tfp - t0p,"] "
