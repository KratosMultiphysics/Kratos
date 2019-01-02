from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

# Import shutil to manage file copying
import shutil

class FracturePropagationUtility(object):

    def __init__(self, model, order_processes_initialization):

        self.model = model
        self.order_processes_initialization = order_processes_initialization

        import platform
        if platform.system()=="Windows":
            self.execute_gid = "gid"
        else:
            self.execute_gid = "./gid"

        # Define FracturesData
        with open("FracturesData.json",'r') as parameter_file:
            self.FracturesData = KratosMultiphysics.Parameters(parameter_file.read())

        # Define control variables
        self.propagation_frequency = self.FracturesData["fracture_data"]["propagation_frequency"].GetInt()
        self.step_count = 0
        self.propagation_count = self.propagation_frequency
        self.remesh_count = 0

        # Define paths
        self.problem_path = os.getcwd()
        self.gid_path = self.FracturesData["fracture_data"]["gid_path"].GetString()
        self.orig_state_path = os.path.join(str(self.problem_path),"OriginalState")
        self.last_state_path = os.path.join(str(self.problem_path),"LastState")

    def Initialize(self, parameters):

        self.domain_size = parameters["solver_settings"]["domain_size"].GetInt()
        if self.domain_size == 2:
            self.PropagationUtility = KratosPoro.FracturePropagation2DUtilities()
            self.tcl_proc = "Poromechanics_Application::PropagateFractures2D"
        else:
            self.PropagationUtility = KratosPoro.FracturePropagation3DUtilities()
            self.tcl_proc = "Poromechanics_Application::PropagateFractures3D"

        self.move_mesh_flag = parameters["solver_settings"]["move_mesh_flag"].GetBool()

        # Create the file containing a list with all post.bin files
        self.problem_name = parameters["problem_data"]["problem_name"].GetString()
        all_list_filename = str(self.problem_name)+"_all.post.lst"
        all_list_file = open(all_list_filename,'w')
        all_list_file.write("Multiple\n")
        all_list_file.close()
        # Save files of the original state
        self.SaveInitialProblemFiles()

        # Modify model_part_name of the input parameters
        self.original_model_part_name = parameters["solver_settings"]["model_part_name"].GetString()
        self.model_part_number = 0
        self.model_part_name = str(self.original_model_part_name) + '_' + str(self.model_part_number)
        parameters["solver_settings"]["model_part_name"].SetString(self.model_part_name)

        if parameters.Has("processes"):
            for name, value in parameters["processes"].items():
                value = self.UpdateModelPartNames(value)

        if parameters.Has("output_processes"):
            for name, value in parameters["output_processes"].items():
                value = self.UpdateModelPartNames(value)

        return parameters

    def UpdateModelPartNames(self,process_parameters):
        for p in range(process_parameters.size()):
            old_name = process_parameters[p]["Parameters"]["model_part_name"].GetString()
            a,b=old_name.split(".")
            new_name = str(self.original_model_part_name) + '_' + str(self.model_part_number) + '.' + str(b)
            process_parameters[p]["Parameters"]["model_part_name"].SetString(new_name)
        return process_parameters

    def SaveInitialProblemFiles(self):

        # Create two folders
        if not os.path.isdir(self.orig_state_path):
            os.mkdir(str(self.orig_state_path))
        else:
            shutil.rmtree(str(self.orig_state_path), ignore_errors=True)
            os.mkdir(str(self.orig_state_path))

        if not os.path.isdir(self.last_state_path):
            os.mkdir(str(self.last_state_path))
        else:
            shutil.rmtree(str(self.last_state_path), ignore_errors=True)
            os.mkdir(str(self.last_state_path))

        # Save list of files names
        self.list_of_files_names = []
        filename = str(self.problem_name)+".geo"
        self.list_of_files_names.append(filename)
        filename = str(self.problem_name)+".lin"
        self.list_of_files_names.append(filename)
        filename = str(self.problem_name)+".prj"
        self.list_of_files_names.append(filename)
        filename = str(self.problem_name)+".tree"
        self.list_of_files_names.append(filename)
        filename = str(self.problem_name)+".msh"
        self.list_of_files_names.append(filename)
        filename = str(self.problem_name)+".cnd"
        self.list_of_files_names.append(filename)
        filename = str(self.problem_name)+".prb"
        self.list_of_files_names.append(filename)
        filename = str(self.problem_name)+".mdpa"
        self.list_of_files_names.append(filename)
        filename = "ProjectParameters.json"
        self.list_of_files_names.append(filename)
        filename = "FracturesData.json"
        self.list_of_files_names.append(filename)

        for filename in self.list_of_files_names:
            filepath = os.path.join(str(self.problem_path),str(filename))
            shutil.copy(str(filepath), str(self.orig_state_path))
            shutil.copy(str(filepath), str(self.last_state_path))

    def IsPropagationStep(self):

        self.step_count += 1
        if self.step_count == self.propagation_count:
            self.propagation_count += self.propagation_frequency
            return True
        else:
            return False

    def CheckPropagation(self,solver,list_of_processes,list_of_output_processes):

        main_model_part = self.model.GetModelPart(self.model_part_name)

        # Check fracture propagation
        propagate_fractures = False
        propagate_fractures = self.PropagationUtility.CheckFracturePropagation(self.FracturesData,main_model_part,self.move_mesh_flag)

        # Generate new fractures if needed
        if propagate_fractures:

            self.remesh_count += 1

            # Overwrite current problem files with last state files
            for filename in self.list_of_files_names:
                filepath = os.path.join(str(self.last_state_path),str(filename))
                shutil.copy(str(filepath), str(self.problem_path))

            # Call GiD to generate new mesh
            import subprocess
            os.chdir(self.gid_path)
            subprocess.call(str(self.execute_gid) + " -n -t \"" + str(self.tcl_proc) + "\" " + str(self.problem_path),shell=True)
            # subprocess.call(str(self.execute_gid) + " -t \"" + str(self.tcl_proc) + "\" " + str(self.problem_path),shell=True)
            os.chdir(self.problem_path)

            # Overwrite last state files with new problem files
            for filename in self.list_of_files_names:
                filepath = os.path.join(str(self.problem_path),str(filename))
                shutil.copy(str(filepath), str(self.last_state_path))

            # Update FracturesData
            with open("FracturesData.json",'r') as parameter_file:
                self.FracturesData = KratosMultiphysics.Parameters(parameter_file.read())

            solver,list_of_processes,list_of_output_processes = self.UpdateSolverAndProcesses(solver,
                                                                                            list_of_processes,
                                                                                            list_of_output_processes)

            # Overwrite current problem files with original state files
            for filename in self.list_of_files_names:
                filepath = os.path.join(str(self.orig_state_path),str(filename))
                shutil.copy(str(filepath), str(self.problem_path))


        return solver,list_of_processes,list_of_output_processes

    def UpdateSolverAndProcesses(self,solver,list_of_processes,list_of_output_processes):

        old_main_model_part = self.model.GetModelPart(self.model_part_name)
        ### Finalize Old Model ---------------------------------------------------------------------------------------

        # Finalizing output files
        list_of_output_processes.ExecuteFinalize()

        for process in list_of_processes:
            process.ExecuteFinalize()

        # Finalizing strategy
        solver.Clear()

        # Save old .post.list file
        all_list_filename = str(self.problem_name)+"_all.post.lst"
        all_list_file = open(all_list_filename,'a')
        partial_list_filename = str(self.problem_name)+".post.lst"
        with open(partial_list_filename) as partial_list_file:
            next(partial_list_file)
            for line in partial_list_file:
                    all_list_file.write(line)
        all_list_file.close()
        # Save old .time file
        original_filename = str(self.problem_name)+".time"
        original_filepath = os.path.join(str(self.problem_path),str(original_filename))
        new_filename = str(self.problem_name)+"_"+str(self.remesh_count)+"info.time"
        new_filepath = os.path.join(str(self.problem_path),str(new_filename))
        shutil.copy(str(original_filepath), str(new_filepath))

        ### Generate New Model ---------------------------------------------------------------------------------------

        # Parsing the parameters
        with open("ProjectParameters.json",'r') as parameter_file:
            ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

        ## Model part ------------------------------------------------------------------------------------------------

        # Defining the model part
        # TODO: I need to change the model_part_name of the processes
        self.model_part_number = self.model_part_number + 1
        new_model_part_name = str(self.original_model_part_name) + '_' + str(self.model_part_number)

        new_model_part = self.model.CreateModelPart(new_model_part_name,2)

        new_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
        new_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, old_main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME])
        new_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, old_main_model_part.ProcessInfo[KratosPoro.TIME_UNIT_CONVERTER])

        # Construct the solver (main setting methods are located in the solver_module)
        solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(new_model_part, ProjectParameters["solver_settings"])

        # Add problem variables
        solver.AddVariables()

        # Read model_part (note: the buffer_size is set here)
        solver.ImportModelPart()

        # Add degrees of freedom
        solver.AddDofs()

        # Print model_part
        echo_level = ProjectParameters["solver_settings"]["echo_level"].GetInt()
        if(echo_level > 1):
            print(new_model_part)

        ## Initialize ------------------------------------------------------------------------------------------------

        # Construct processes to be applied
        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
        list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )

        # Initialize processes
        for process in list_of_processes:
            process.ExecuteInitialize()

        # Set TIME
        new_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, old_main_model_part.ProcessInfo[KratosMultiphysics.TIME])

        # Initialize GiD I/O
        computing_model_part = solver.GetComputingModelPart()
        output_settings = ProjectParameters["output_configuration"]
        from gid_output_process import GiDOutputProcess
        list_of_output_processes = GiDOutputProcess(computing_model_part,self.problem_name,output_settings)
        list_of_output_processes.ExecuteInitialize()

        # Initialize the solver
        solver.Initialize()

        # Initialize the strategy before the mapping
        solver.InitializeStrategy()

        # ExecuteBeforeSolutionLoop
        for process in list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Set results when they are written in a single file (only multiplefiles for the moment)
        list_of_output_processes.ExecuteBeforeSolutionLoop()

        ### Mapping between old and new model parts ------------------------------------------------------------------

        self.PropagationUtility.MappingModelParts(self.FracturesData,old_main_model_part,new_model_part,self.move_mesh_flag)

        # set ARC_LENGTH_LAMBDA and ARC_LENGTH_RADIUS_FACTOR and update loads
        if ProjectParameters["solver_settings"]["strategy_type"].GetString() == "arc_length":
            new_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_LAMBDA, old_main_model_part.ProcessInfo[KratosPoro.ARC_LENGTH_LAMBDA])
            new_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_RADIUS_FACTOR, old_main_model_part.ProcessInfo[KratosPoro.ARC_LENGTH_RADIUS_FACTOR])
            solver._UpdateLoads()

        # delete old model_part
        old_model_part_number = self.model_part_number - 1
        old_model_part_name = str(self.original_model_part_name) + '_' + str(old_model_part_number)
        self.model.DeleteModelPart(old_model_part_name)

        # Check new mesh
        #IsConverged = solver._CheckConvergence()

        return solver,list_of_processes,list_of_output_processes

    def Finalize(self):

        # Save old .post.list file
        all_list_filename = str(self.problem_name)+"_all.post.lst"
        all_list_file = open(all_list_filename,'a')
        partial_list_filename = str(self.problem_name)+".post.lst"
        with open(partial_list_filename) as partial_list_file:
            next(partial_list_file)
            for line in partial_list_file:
                    all_list_file.write(line)
        all_list_file.close()
