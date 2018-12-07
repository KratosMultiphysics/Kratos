from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.CompressiblePotentialFlowApplication import *
import potential_flow_solver
import process_factory

class PotentialFlowAnalysisTest:

    def __init__(self, model, project_parameters):
        """The constructor of the PotentialAnalysis-Object.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        project_parameters -- The self.project_parameters used
        """
        if (type(model) != KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(project_parameters) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

        ## Fluid model part definition
        self.main_model_part = model.CreateModelPart(project_parameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, project_parameters["problem_data"]["domain_size"].GetInt())

        ###TODO replace this "model" for real one once available
        self.model = {self.project_parameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        ## Solver construction
        self.solver = potential_flow_solver.CreateSolver(self.main_model_part, project_parameters["solver_settings"])

        self.solver.AddVariables()

    def Initialize(self):
        ## Read the model - note that SetBufferSize is done here
        self.solver.ImportModelPart()

        ## Add AddDofs
        self.solver.AddDofs()

        ## Initialize GiD  I/O
        if (self.parallel_type == "OpenMP"):
            from gid_output_process import GiDOutputProcess
            self.gid_output = GiDOutputProcess(self.solver.GetComputingModelPart(),
                                          self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                          self.project_parameters["output_configuration"])

        self.gid_output.ExecuteInitialize()

        ##TODO: replace MODEL for the Kratos one ASAP
        ## Get the list of the skin submodel parts in the object Model
        for i in range(self.project_parameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = self.project_parameters["solver_settings"]["skin_parts"][i].GetString()
            self.model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})

        ## Get the list of the no-skin submodel parts in the object Model (results processes and no-skin conditions)
        for i in range(self.project_parameters["solver_settings"]["no_skin_parts"].size()):
            no_skin_part_name = self.project_parameters["solver_settings"]["no_skin_parts"][i].GetString()
            self.model.update({no_skin_part_name: self.main_model_part.GetSubModelPart(no_skin_part_name)})


        ## Print model_part and properties
        if(self.echo_level > 1):
            print("")
            print(self.main_model_part)
            for properties in self.main_model_part.Properties:
                print(properties)

        ## Processes construction
        # "list_of_processes" contains all the processes already constructed (boundary conditions and initial conditions)
        self.list_of_processes = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.project_parameters["boundary_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.project_parameters["auxiliar_process_list"] )

        if(self.echo_level > 1):
            for process in self.list_of_processes:
                print(process)

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        ## Solver initialization
        self.solver.Initialize()

        #TODO: think if there is a better way to do this
        self.fluid_model_part = self.solver.GetComputingModelPart()

        self.gid_output.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Writing the full self.project_parameters file before solving
        if ((self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0)) and (self.echo_level > 0):
            f = open("self.project_parametersOutput.json", 'w')
            f.write(self.project_parameters.PrettyPrintJsonString())
            f.close()

    def RunSolution(self):
        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        self.gid_output.ExecuteInitializeSolutionStep()

        self.solver.Solve()

        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        self.gid_output.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

        if self.gid_output.IsOutputStep():
            self.gid_output.PrintOutput()

        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()

        for process in self.list_of_processes:
            process.ExecuteFinalize()

        self.gid_output.ExecuteFinalize()

    def Run(self):
        self.Initialize()
        self.RunSolution()
