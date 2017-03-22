from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMeshing
import KratosMultiphysics.ExternalSolversApplication as KratosExternalSolvers
import KratosMultiphysics.ALEApplication as KratosALE

import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteMeshTest(KratosUnittest.TestCase):

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        ## Solver construction
        mesh_solver_type = self.ProjectParameters["solver_settings"]["solver_type"].GetString()

        # Creation an empty parameters object to set the solver defaults
        mesh_solver_parameters = KratosMultiphysics.Parameters("""{}""")

        # Mesh solver construction
        if(mesh_solver_type == "Laplacian"):
          import mesh_solver_laplacian as mesh_solver
          self.mesh_sol = mesh_solver.MeshSolverLaplacian(self.main_model_part, mesh_solver_parameters)

        elif(mesh_solver_type == "StructuralSimilarity"):
          import mesh_solver_structural_similarity as mesh_solver
          self.mesh_sol = mesh_solver.MeshSolverStructuralSimilarity(self.main_model_part, mesh_solver_parameters)

        elif(mesh_solver_type == "Ballvertex"):
          import mesh_solver_ballvertex as mesh_solver
          self.mesh_sol = mesh_solver.MeshSolverBallvertex(self.main_model_part, mesh_solver_parameters)

        else:
          raise NameError("Selected solver has not been implemented!")

        ## Mesh solver variables addition
        self.mesh_sol.AddVariables()

        ## Read the model - note that SetBufferSize is done here
        input_file_name = self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()   # introducing input file name
        KratosMultiphysics.ModelPartIO(input_file_name).ReadModelPart(self.main_model_part)                                  # reading the fluid part
        self.main_model_part.SetBufferSize(3)                                                                                # setting up the buffer size

        ## Mesh solver DOFs addition
        self.mesh_sol.AddDofs()

        ## Mesh solver initialize
        self.mesh_sol.Initialize()

        ## Initialize GiD  I/O
        self.output_flag = False
        if (self.output_flag == True):
            from gid_output_process import GiDOutputProcess
            self.gid_output = GiDOutputProcess(self.mesh_sol.GetComputingModelPart(),
                                               self.ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                               self.ProjectParameters["output_configuration"])

            self.gid_output.ExecuteInitialize()

        ## Get the fluid submodel part in the object Model
        fluid_part_name = ProjectParameters["solver_settings"]["volume_model_part_name"].GetString()
        Model.update({fluid_part_name: self.main_model_part.GetSubModelPart(fluid_part_name)})
        ## Get the list of the skin submodel parts in the object Model
        for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
            Model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})

        ## Processes construction
        import process_factory
        self.list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()


    def Solve(self):

        ## Stepping and time settings
        Dt = self.ProjectParameters["problem_data"]["delta_time"].GetDouble()
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = 0.0
        step = 0

        if (self.output_flag == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        while(time <= end_time):

            time = time + Dt
            step = step + 1
            self.main_model_part.CloneTimeStep(time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()

            if (self.output_flag == True):
                self.gid_output.ExecuteInitializeSolutionStep()

            self.mesh_sol.Solve()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()

            if (self.output_flag == True):
                self.gid_output.ExecuteFinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            if (self.output_flag == True):
                self.gid_output.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

        for process in self.list_of_processes:
            process.ExecuteFinalize()
