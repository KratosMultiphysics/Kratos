from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Import system python
import os
# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication  as KratosSolid
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
# Import shutil to manage file copying
import shutil

class FracturePropagationUtility:

    def __init__(self,domain_size,problem_name,move_mesh_flag):

        # Construct the utility
        self.domain_size = domain_size
        if domain_size==2:
            self.PropagationUtility = KratosPoro.FracturePropagation2DUtilities()
            self.tcl_proc = "Poromechanics_Application::PropagateFractures2D"
        else:
            self.PropagationUtility = KratosPoro.FracturePropagation3DUtilities()
            self.tcl_proc = "Poromechanics_Application::PropagateFractures3D"

        import platform
        if platform.system()=="Windows":
            self.execute_gid = "gid"
        else:
            self.execute_gid = "./gid"

        self.move_mesh_flag = move_mesh_flag

        # Define FracturesData
        with open("FracturesData.json",'r') as parameter_file:
            self.FracturesData = KratosMultiphysics.Parameters(parameter_file.read())

        # Define control variables
        self.propagation_frequency = self.FracturesData["fracture_data"]["propagation_frequency"].GetInt()
        self.step_count = 0
        self.propagation_count = self.propagation_frequency
        self.remesh_count = 0

        # Define names and paths
        self.problem_name = problem_name
        self.problem_path = os.getcwd()
        self.gid_path = self.FracturesData["fracture_data"]["gid_path"].GetString()
        self.orig_state_path = os.path.join(str(self.problem_path),"OriginalState")
        self.last_state_path = os.path.join(str(self.problem_path),"LastState")

        # Create the file containing a list with all post.bin files
        all_list_filename = str(self.problem_name)+"_all.post.lst"
        all_list_file = open(all_list_filename,'w')
        all_list_file.write("Multiple\n")
        all_list_file.close()

        # Save files of the original state
        self.SaveInitialProblemFiles()

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

    def CheckPropagation(self,main_model_part,solver,list_of_processes,gid_output):

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

            main_model_part,solver,list_of_processes,gid_output = self.GenereateNewModelPart(main_model_part,
                                                                                             solver,
                                                                                             list_of_processes,
                                                                                             gid_output)

            # Overwrite current problem files with original state files
            for filename in self.list_of_files_names:
                filepath = os.path.join(str(self.orig_state_path),str(filename))
                shutil.copy(str(filepath), str(self.problem_path))


        return main_model_part,solver,list_of_processes,gid_output

    def GenereateNewModelPart(self,main_model_part,solver,list_of_processes,gid_output):

        ### Finalize Old Model ---------------------------------------------------------------------------------------

        # Finalizing output files
        gid_output.ExecuteFinalize()

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

        # Save previous model_part
        main_model_part_old = main_model_part

        ### Generate New Model ---------------------------------------------------------------------------------------

        # Parsing the parameters
        with open("ProjectParameters.json",'r') as parameter_file:
            ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

        ## Model part ------------------------------------------------------------------------------------------------

        # Defining the model part
        main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["solver_settings"]["model_part_name"].GetString())
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, main_model_part_old.ProcessInfo[KratosMultiphysics.DELTA_TIME])
        main_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, main_model_part_old.ProcessInfo[KratosPoro.TIME_UNIT_CONVERTER])

        # Construct the solver (main setting methods are located in the solver_module)
        solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

        # Add problem variables
        solver.AddVariables()

        # Read model_part (note: the buffer_size is set here)
        solver.ImportModelPart()

        # Add degrees of freedom
        solver.AddDofs()

        # Creation of Kratos model
        PoroModel = KratosMultiphysics.Model()
        PoroModel.AddModelPart(main_model_part)

        # Build sub_model_parts (save the list of the submodel part in the object Model)
        for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            PoroModel.AddModelPart(main_model_part.GetSubModelPart(part_name))

        # Print model_part
        echo_level = ProjectParameters["solver_settings"]["echo_level"].GetInt()
        if(echo_level > 1):
            print(main_model_part)

        ## Initialize ------------------------------------------------------------------------------------------------

        # Construct processes to be applied
        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(PoroModel).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
        list_of_processes += process_factory.KratosProcessFactory(PoroModel).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )

        # Initialize processes
        for process in list_of_processes:
            process.ExecuteInitialize()

        # Set TIME
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, main_model_part_old.ProcessInfo[KratosMultiphysics.TIME])

        # Initialize GiD I/O
        computing_model_part = solver.GetComputingModelPart()
        output_settings = ProjectParameters["output_configuration"]
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(computing_model_part,self.problem_name,output_settings)
        gid_output.ExecuteInitialize()

        # Initialize the solver
        solver.Initialize()

        # Initialize the strategy before the mapping
        solver.InitializeStrategy()

        # ExecuteBeforeSolutionLoop
        for process in list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Set results when they are written in a single file (only multiplefiles for the moment)
        gid_output.ExecuteBeforeSolutionLoop()

        ### Mapping between old and new model parts ------------------------------------------------------------------

        self.PropagationUtility.MappingModelParts(self.FracturesData,main_model_part_old,main_model_part,self.move_mesh_flag)

        # set ARC_LENGTH_LAMBDA and ARC_LENGTH_RADIUS_FACTOR and update loads
        if ProjectParameters["solver_settings"]["strategy_type"].GetString() == "arc_length":
            main_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_LAMBDA, main_model_part_old.ProcessInfo[KratosPoro.ARC_LENGTH_LAMBDA])
            main_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_RADIUS_FACTOR, main_model_part_old.ProcessInfo[KratosPoro.ARC_LENGTH_RADIUS_FACTOR])
            solver._UpdateLoads()

        # delete auxiliary model_part
        del main_model_part_old

        # Check new mesh
        #IsConverged = solver._CheckConvergence()

        return main_model_part,solver,list_of_processes,gid_output

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
