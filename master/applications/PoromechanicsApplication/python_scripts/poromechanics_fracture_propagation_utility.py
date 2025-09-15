
# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

# Import shutil to manage file copying
import shutil

from importlib import import_module

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
        self.model_part_number = 0

        # Define paths
        self.problem_path = os.getcwd()
        self.gid_path = self.FracturesData["fracture_data"]["gid_path"].GetString()
        self.orig_state_path = os.path.join(str(self.problem_path),"OriginalState")
        self.last_state_path = os.path.join(str(self.problem_path),"LastState")

    def Initialize(self, parameters):

        domain_size = parameters["solver_settings"]["domain_size"].GetInt()
        if domain_size == 2:
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

        # Update parameters with the new model part name
        self.original_model_part_name = parameters["solver_settings"]["model_part_name"].GetString()
        model_part_name = str(self.original_model_part_name) + '_' + str(self.model_part_number)

        parameters["solver_settings"]["model_part_name"].SetString(model_part_name)

        if parameters.Has("processes"):
            for name, value in parameters["processes"].items():
                value = self.UpdateModelPartNames(value)

        if parameters.Has("output_processes"):
            for name, value in parameters["output_processes"].items():
                value = self.UpdateModelPartNames(value)

        # Overwrite materials with the new model part name
        with open("PoroMaterials.json",'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())
        for p in range(materials["properties"].size()):
            old_name = materials["properties"][p]["model_part_name"].GetString()
            a,b=old_name.split(".")
            new_name = str(self.original_model_part_name) + '_' + str(self.model_part_number) + '.' + str(b)
            materials["properties"][p]["model_part_name"].SetString(new_name)
        with open("PoroMaterials.json", 'w') as outfile:
            outfile.write(materials.PrettyPrintJsonString())

        return parameters

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
        filename = "PoroMaterials.json"
        self.list_of_files_names.append(filename)
        filename = "FracturesData.json"
        self.list_of_files_names.append(filename)

        for filename in self.list_of_files_names:
            filepath = os.path.join(str(self.problem_path),str(filename))
            shutil.copy(str(filepath), str(self.orig_state_path))
            shutil.copy(str(filepath), str(self.last_state_path))

    def UpdateModelPartNames(self,process_parameters):
        for p in range(process_parameters.size()):
            old_name = process_parameters[p]["Parameters"]["model_part_name"].GetString()
            a,b=old_name.split(".")
            new_name = str(self.original_model_part_name) + '_' + str(self.model_part_number) + '.' + str(b)
            process_parameters[p]["Parameters"]["model_part_name"].SetString(new_name)
        return process_parameters

    def IsPropagationStep(self):

        self.step_count += 1
        if self.step_count == self.propagation_count:
            self.propagation_count += self.propagation_frequency
            return True
        else:
            return False

    def CheckPropagation(self,solver,list_of_processes,list_of_output_processes):

        model_part_name = str(self.original_model_part_name) + '_' + str(self.model_part_number)
        main_model_part = self.model.GetModelPart(model_part_name)

        # Check fracture propagation
        propagate_fractures = False
        propagate_fractures = self.PropagationUtility.CheckFracturePropagation(self.FracturesData,main_model_part,self.move_mesh_flag)

        # Generate new fractures if needed
        if propagate_fractures:

            self.model_part_number += 1

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

        ## Finalize old solver and processes

        # Finalize list_of_processes and list_of_output_processes (included in list_of_processes)
        for process in list_of_processes:
            process.ExecuteFinalize()

        solver.Finalize()

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

        ## Create new solver and processes

        # Parsing the parameters
        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # Update parameters with the new model part name
        new_model_part_name = str(self.original_model_part_name) + '_' + str(self.model_part_number)
        parameters["solver_settings"]["model_part_name"].SetString(new_model_part_name)

        if parameters.Has("processes"):
            for name, value in parameters["processes"].items():
                value = self.UpdateModelPartNames(value)

        if parameters.Has("output_processes"):
            for name, value in parameters["output_processes"].items():
                value = self.UpdateModelPartNames(value)

        # Overwrite materials with the new model part name
        with open("PoroMaterials.json",'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())
        for p in range(materials["properties"].size()):
            old_name = materials["properties"][p]["model_part_name"].GetString()
            a,b=old_name.split(".")
            new_name = str(self.original_model_part_name) + '_' + str(self.model_part_number) + '.' + str(b)
            materials["properties"][p]["model_part_name"].SetString(new_name)
        with open("PoroMaterials.json", 'w') as outfile:
            outfile.write(materials.PrettyPrintJsonString())

        # Create new solver (and new_model_part)
        python_module_name = "KratosMultiphysics.PoromechanicsApplication"
        full_module_name = python_module_name + "." + parameters["solver_settings"]["solver_type"].GetString()
        solver_module = import_module(full_module_name)
        new_solver = solver_module.CreateSolver(self.model, parameters["solver_settings"])

        new_solver.AddVariables()
        new_solver.ImportModelPart()
        new_solver.PrepareModelPart()
        new_solver.AddDofs()

        # Create new_list_of_processes
        new_list_of_processes = []
        new_list_of_output_processes = []

        from KratosMultiphysics.process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)

        if parameters.Has("processes"):
            processes_params = parameters["processes"]

            # first initialize the processes that depend on the order
            for processes_names in self.order_processes_initialization:
                if processes_params.Has(processes_names):
                    new_list_of_processes += factory.ConstructListOfProcesses(processes_params[processes_names])

            # then initialize the processes that don't depend on the order
            for name, value in processes_params.items():
                if not name in self.order_processes_initialization:
                    new_list_of_processes += factory.ConstructListOfProcesses(value)

        order_output_processes_initialization = []
        if parameters.Has("output_processes"):
            processes_params = parameters["output_processes"]

            # first initialize the processes that depend on the order
            for processes_names in order_output_processes_initialization:
                if processes_params.Has(processes_names):
                    new_list_of_output_processes += factory.ConstructListOfProcesses(processes_params[processes_names])

            # then initialize the processes that don't depend on the order
            for name, value in processes_params.items():
                if not name in order_output_processes_initialization:
                    new_list_of_output_processes += factory.ConstructListOfProcesses(value)

        new_list_of_processes.extend(new_list_of_output_processes)

        # Initialize processes and solver
        for process in new_list_of_processes:
            process.ExecuteInitialize()

        new_solver.Initialize()

        new_solver.Check()
        for process in new_list_of_processes:
            process.Check()

        for process in new_list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Mapping between old and new model parts

        old_model_part_name = str(self.original_model_part_name) + '_' + str(self.model_part_number-1)

        old_main_model_part = self.model.GetModelPart(old_model_part_name)
        new_model_part = self.model.GetModelPart(new_model_part_name)

        new_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, old_main_model_part.ProcessInfo[KratosMultiphysics.TIME])
        new_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, old_main_model_part.ProcessInfo[KratosMultiphysics.STEP])

        self.PropagationUtility.MappingModelParts(self.FracturesData,old_main_model_part,new_model_part,self.move_mesh_flag)

        # Set ARC_LENGTH_LAMBDA and ARC_LENGTH_RADIUS_FACTOR and update loads
        if parameters["solver_settings"]["strategy_type"].GetString() == "arc_length":
            new_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_LAMBDA, old_main_model_part.ProcessInfo[KratosPoro.ARC_LENGTH_LAMBDA])
            new_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_RADIUS_FACTOR, old_main_model_part.ProcessInfo[KratosPoro.ARC_LENGTH_RADIUS_FACTOR])
            new_solver._UpdateLoads()

        ## Delete old model_part and replace solver and processes

        self.model.DeleteModelPart(old_model_part_name)
        solver = new_solver
        list_of_processes = new_list_of_processes
        list_of_output_processes = new_list_of_output_processes

        # Check new mesh
        # IsConverged = new_solver._CheckConvergence()

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
