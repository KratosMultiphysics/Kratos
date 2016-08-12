##################################################################################################
# GENERAL MAIN OF THE ANALISYS
##################################################################################################

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# time control starts
from time import *
print(ctime())

# import additional libraries
import math
import csv
import os
from itertools import islice

# including kratos path
from KratosMultiphysics import *

# including Applications paths
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

# import the python utilities:
import restart_utility as restart_utils
#import gid_output_utility as gid_utils
from gid_output import GiDOutput
import conditions_python_utility as condition_utils
import list_files_python_utility as files_utils
import time_operation_utility as operation_utils

##################################################################################################

class Toolbox:

    # ********************************************************************************************
    def __init__(self, parameters):

        # Initialize model part
        self.model_part = "not_read_yet"

        # Store parameters
        self.parameters = parameters

        # defining the number of threads:
        num_threads = self.parameters.NumberofThreads
        self.SetParallelSize(num_threads)
        
        # Container to store boundary deformation
        self.boundary_deformation = {}

        # Read CSM mesh
        self.read_mesh()

    # ********************************************************************************************
    def read_mesh(self):

        print("\n> Start reading CSM mesh...")

        # Take time
        t0p = clock()

        # defining the type, the name and the path of the problem:
        problem_type = self.parameters.ProblemType
        self.problem_name = self.parameters.problem_name
        # Get current working directory
        current_directory = os.getcwd()
        problem_path = current_directory

        # defining a model part
        self.model_part = ModelPart("SolidDomain")

        # defining the model size to scale
        length_scale = 1.0

        # --DEFINE MAIN SOLVER START

        SolverSettings = self.parameters.SolverSettings

        # import solver file
        solver_constructor = __import__(SolverSettings.solver_type)

        # construct the solver
        self.main_step_solver = solver_constructor.CreateSolver(self.model_part, SolverSettings)

        # --DEFINE MAIN SOLVER END


        # --READ AND SET MODEL FILES

        # set the restart of the problem
        restart_step = self.parameters.Restart_Step
        self.problem_restart = restart_utils.RestartUtility(self.model_part, problem_path, self.problem_name)

        # set the results file list of the problem (managed by the problem_restart and gid_print)
        print_lists = self.parameters.PrintLists
        output_mode = self.parameters.GidOutputConfiguration.GiDPostMode
        self.list_files = files_utils.ListFilesUtility(problem_path, self.problem_name, print_lists, output_mode)
        self.list_files.Initialize(self.parameters.file_list)

        # --READ AND SET MODEL FILES END


        # --DEFINE CONDITIONS START
        incr_disp = self.parameters.Incremental_Displacement
        incr_load = self.parameters.Incremental_Load
        rotation_dofs = SolverSettings.RotationDofs
        self.conditions = condition_utils.ConditionsUtility(self.model_part, self.parameters.domain_size, incr_disp, incr_load, rotation_dofs)

        # --DEFINE CONDITIONS END

        # --CONFIGURATIONS END
        # ----------------------------------------------------------------#

        # --START SOLUTION
        #
        # initialize problem : load restart or initial start
        self.load_restart = self.parameters.LoadRestart
        self.save_restart = self.parameters.SaveRestart

        # set buffer size
        self.buffer_size = 3

        # define problem variables:
        solver_constructor.AddVariables(self.model_part, SolverSettings)

        # --- READ MODEL
        if(self.load_restart == False):

            # remove results, restart, graph and list previous files
            self.problem_restart.CleanPreviousFiles()
            self.list_files.RemoveListFiles()

            # reading the model
            self.model_part_io = ModelPartIO(self.problem_name)
            self.model_part_io.ReadModelPart(self.model_part)

            # set the buffer size
            self.model_part.SetBufferSize(self.buffer_size)
            # Note: the buffer size should be set once the mesh is read for the first time

            # set the degrees of freedom
            solver_constructor.AddDofs(self.model_part, SolverSettings)

            # set the constitutive law
            import constitutive_law_python_utility as constitutive_law_utils

            constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.model_part, self.parameters.domain_size);
            constitutive_law.Initialize();

        else:

            # reading the model from the restart file
            self.problem_restart.Load(restart_step);

            # remove results, restart, graph and list posterior files
            self.problem_restart.CleanPosteriorFiles(restart_step)
            self.list_files.ReBuildListFiles()

        # initialize mesh motion solver
        self.initialize_solver()

        # measure process time
        tfp = clock()

        # --INITIALIZE GID I/0
        self.gid_io = GiDOutput(self.problem_name,
                           self.parameters.GidOutputConfiguration.VolumeOutput,
                           self.parameters.GidOutputConfiguration.GiDPostMode,
                           self.parameters.GidOutputConfiguration.GiDMultiFileFlag,
                           self.parameters.GidOutputConfiguration.GiDWriteMeshFlag,
                           self.parameters.GidOutputConfiguration.GiDWriteConditionsFlag)

        # --GID OUTPUT OPTIONS END

        print("> Reading CSM mesh completed [Reading time = ", tfp - t0p, "]\n")        

    # ********************************************************************************************
    def initialize_solver(self):

        # set delta time in process info
        self.model_part.ProcessInfo[DELTA_TIME] = self.parameters.time_step

        # solver initialize
        self.main_step_solver.Initialize()
        self.main_step_solver.SetRestart(self.load_restart) #calls strategy initialize if no restart

        # define time steps and loop range of steps
        if(self.load_restart):
            self.buffer_size = 0
        else:
            self.model_part.ProcessInfo[TIME]                = 0
            self.model_part.ProcessInfo[TIME_STEPS]          = 0
            self.conditions.Initialize(self.model_part.ProcessInfo[DELTA_TIME]);

    # ********************************************************************************************
    def initialize_output(self):

        # writing a single file
        self.gid_io.initialize_results(self.model_part)       

        # writing a initial state results file
        self.gid_io.write_results(0,self.model_part,self.parameters.nodal_results,self.parameters.gauss_points_results)
   
    # ********************************************************************************************
    def finalize_output(self):

        # writing a single file
        self.gid_io.finalize_results()

    # ********************************************************************************************
    def solve(self,nodal_loads,itr):

        print("\n> Start solving CSM...")

        # --TIME INTEGRATION

        # Take time
        t0p = clock()

        # initialize step operations
        self.starting_step  = 0
        self.starting_time  = 0
        self.ending_time    = self.parameters.end_time

        #initialize time integration variables
        self.current_time = self.starting_time
        self.current_step = self.starting_step

        # filling the buffer
        for step in range(0,self.buffer_size):

          self.model_part.CloneTimeStep(self.current_time)
          self.model_part.ProcessInfo[TIME_STEPS] = step-self.buffer_size

        # store previous time step
        self.model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = self.model_part.ProcessInfo[DELTA_TIME]

        # Advance
        self.current_time = self.current_time + self.model_part.ProcessInfo[DELTA_TIME]
        self.current_step = self.current_step + 1
        self.model_part.CloneTimeStep(self.current_time)
        self.model_part.ProcessInfo[TIME] = self.current_time

        # Apply given load boundary condition
        self.apply_load_condition(nodal_loads)

        # Solve mesh
        self.main_step_solver.Solve()

        # incremental load
        self.conditions.SetIncrementalLoad(self.current_step, self.model_part.ProcessInfo[DELTA_TIME]);

        # print the results at the end of the step
        self.gid_io.write_results(itr,self.model_part,self.parameters.nodal_results,self.parameters.gauss_points_results)

        # self.conditions.RestartImposedDisp()

        # --END

        # measure process time
        tfp = clock()

        print("> Finished solving CSM  [Analysis Time = ", tfp - t0p, "]\n")

        # Return computed displacements
        displacements = {}
        for node in self.model_part.Nodes:
            disp_x = node.GetSolutionStepValue(DISPLACEMENT_X)
            disp_y = node.GetSolutionStepValue(DISPLACEMENT_Y)
            disp_z = node.GetSolutionStepValue(DISPLACEMENT_Z)
            displacements[node.Id] = [disp_x,disp_y,disp_z]
        return displacements

    # ********************************************************************************************   
    def update_mesh(self,X):

        print("> Start writing new CSM mdpa file...")

        # Take time
        t0p = clock()

        # update previous mesh in memory
        for node_ID in X:
            self.model_part.Nodes[node_ID].X = self.model_part.Nodes[node_ID].X0 + X[node_ID][0]
            self.model_part.Nodes[node_ID].Y = self.model_part.Nodes[node_ID].Y0 + X[node_ID][1]
            self.model_part.Nodes[node_ID].Z = self.model_part.Nodes[node_ID].Z0 + X[node_ID][2]
       
        # First we write the new mesh in a temporary file 
        f_tb_created = open(self.problem_name+"_temp.mdpa", 'w')

        # open previous mesh file to read
        f_tb_read = open(self.problem_name+".mdpa",'r')

        # variable to check whether node-section starts
        is_node_section = False

        # scan file until end of file
        keepon = True
        while keepon:

            # read line
            n_lines = 1
            fileslice = islice(f_tb_read,n_lines)
            line = list(fileslice)
            
            # stop if line is empty
            if not line: 
                keepon = False
                break

            # Get actual line content
            line = line[0]

            # Write line if we are not in node section
            if(is_node_section==False): 
                f_tb_created.write(str(line))

            # Check node section starts or ends
            if "Begin Nodes" in line: 
                for node in self.model_part.Nodes:
                    f_tb_created.write(str(node.Id))
                    f_tb_created.write("\t")  
                    f_tb_created.write(str("%.15f"%(node.X)))
                    f_tb_created.write("\t")  
                    f_tb_created.write(str("%.15f"%(node.Y)))
                    f_tb_created.write("\t")
                    f_tb_created.write(str("%.15f"%(node.Z)))
                    f_tb_created.write("\n") 
                is_node_section = True
            elif "End Nodes" in line:
                f_tb_created.write(str(line))
                is_node_section = False
                
        # Close files
        f_tb_read.close()
        f_tb_created.close()

        # Overwrite previous mesh-file with new one
        os.system("cp "+self.parameters.problem_name+"_temp.mdpa "+self.parameters.problem_name+".mdpa")
        os.system("rm -rf "+self.parameters.problem_name+"_temp.mdpa")
         
        # measure process time
        tfp = clock()

        print("> Finished writing mdpa-file!  [Analysis Time = ", tfp - t0p, "]\n")   

        # Re-read CSM mesh after new mesh file was generated
        self.read_mesh()  

    # ********************************************************************************************
    def SetParallelSize(self,num_threads):
        parallel = OpenMPUtils()
        print("Num Threads = ", num_threads)
        parallel.SetNumThreads(int(num_threads))

    # ********************************************************************************************
    def apply_load_condition(self,nodal_loads):

        # Set boundary deformation update
        for node in self.model_part.Nodes:
            load = Vector(3)
            load[0] = nodal_loads[node.Id][0]
            load[1] = nodal_loads[node.Id][1]
            load[2] = nodal_loads[node.Id][2]
            node.SetSolutionStepValue(POINT_LOAD,0,load)

    # ********************************************************************************************

##################################################################################################