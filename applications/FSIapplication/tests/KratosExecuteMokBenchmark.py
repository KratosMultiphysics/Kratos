from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteMokBenchmark(KratosUnittest.TestCase):

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        ## Fluid-Structure model parts definition
        self.structure_main_model_part = ModelPart(self.ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString())
        self.structure_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["structure_solver_settings"]["problem_data"]["domain_size"].GetInt())

        self.fluid_main_model_part = ModelPart(self.ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString())
        self.fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["fluid_solver_settings"]["problem_data"]["domain_size"].GetInt())

        FluidModel = {self.ProjectParameters["fluid_solver_settings"]["problem_data"]["model_part_name"].GetString() : self.fluid_main_model_part}
        SolidModel = {self.ProjectParameters["structure_solver_settings"]["problem_data"]["model_part_name"].GetString() : self.structure_main_model_part}

        ## Solver construction
        solver_module = __import__("partitioned_fsi_solver") # Currently there is only one FSI solver up to date
        self.solver = solver_module.CreateSolver(self.structure_main_model_part, self.fluid_main_model_part, self.ProjectParameters)

        self.solver.AddVariables()

        ## Read the model - note that SetBufferSize is done here
        self.solver.ImportModelPart()

        ## Add AddDofs
        self.solver.AddDofs()

        ## Initialize GiD  I/O
        self.output_post = False # Set this variable to True if it is need to print the results for debugging purposes
        if (self.output_post == True):
            from gid_output_process import GiDOutputProcess

            self.gid_output_structure = GiDOutputProcess(self.solver.structure_solver.GetComputingModelPart(),
                                                         self.ProjectParameters["structure_solver_settings"]["problem_data"]["problem_name"].GetString()+"_structure",
                                                         self.ProjectParameters["structure_solver_settings"]["output_configuration"])

            self.gid_output_fluid = GiDOutputProcess(self.solver.fluid_solver.GetComputingModelPart(),
                                                     self.ProjectParameters["fluid_solver_settings"]["problem_data"]["problem_name"].GetString()+"_fluid",
                                                     self.ProjectParameters["fluid_solver_settings"]["output_configuration"])

            self.gid_output_structure.ExecuteInitialize()
            self.gid_output_fluid.ExecuteInitialize()

        ## Get the list of the skin submodel parts in the object Model (FLUID)
        for i in range(self.ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"].size()):
            skin_part_name = self.ProjectParameters["fluid_solver_settings"]["solver_settings"]["skin_parts"][i].GetString()
            FluidModel.update({skin_part_name: self.fluid_main_model_part.GetSubModelPart(skin_part_name)})

        ## Get the list of the no-skin submodel parts in the object Model (FLUID)
        for i in range(self.ProjectParameters["fluid_solver_settings"]["solver_settings"]["no_skin_parts"].size()):
            no_skin_part_name = self.ProjectParameters["fluid_solver_settings"]["solver_settings"]["no_skin_parts"][i].GetString()
            FluidModel.update({no_skin_part_name: self.fluid_main_model_part.GetSubModelPart(no_skin_part_name)})

        ## Get the list of the initial conditions submodel parts in the object Model (FLUID)
        for i in range(self.ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"].size()):
            initial_cond_part_name = self.ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
            FluidModel.update({initial_cond_part_name: self.fluid_main_model_part.GetSubModelPart(initial_cond_part_name)})

        ## Get the gravity submodel part in the object Model (FLUID)
        for i in range(self.ProjectParameters["fluid_solver_settings"]["gravity"].size()):
            gravity_part_name = self.ProjectParameters["fluid_solver_settings"]["gravity"][i]["Parameters"]["model_part_name"].GetString()
            FluidModel.update({gravity_part_name: self.fluid_main_model_part.GetSubModelPart(gravity_part_name)})

        ## Get the list of the submodel part in the object Model (STRUCTURE)
        for i in range(self.ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = self.ProjectParameters["structure_solver_settings"]["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            SolidModel.update({part_name: self.structure_main_model_part.GetSubModelPart(part_name)})

        ## Processes construction
        import process_factory
        # "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
        # Note that the conditions are firstly constructed. Otherwise, they may overwrite the BCs information.

        # FLUID DOMAIN PROCESSES
        self.list_of_processes = process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( self.ProjectParameters["fluid_solver_settings"]["initial_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( self.ProjectParameters["fluid_solver_settings"]["boundary_conditions_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( self.ProjectParameters["fluid_solver_settings"]["gravity"] )

        # SOLID DOMAIN PROCESSES
        self.list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( self.ProjectParameters["structure_solver_settings"]["constraints_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(SolidModel).ConstructListOfProcesses( self.ProjectParameters["structure_solver_settings"]["loads_process_list"] )

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        # Solver initialization moved after the processes initialization, otherwise the flag INTERFACE is not set
        self.solver.Initialize()

    def Solve(self):

        ## Time settings
        end_time = self.ProjectParameters["fluid_solver_settings"]["problem_data"]["end_time"].GetDouble()
        time = 0.0
        step = 0
        out = 0.0

        if (self.output_post == True):
            self.gid_output_structure.ExecuteBeforeSolutionLoop()
            self.gid_output_fluid.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        while(time <= end_time):

            Dt = (self.solver).ComputeDeltaTime()
            time = time + Dt
            step = step + 1

            # Custom velocity profile for Mok benchmark
            # Note that this is not the original time variation. A linear one has been set to avoid to import Numpy.
            if time <=10.0:
                v_bar = 0.06067*(time/10.0)
            elif time>10.0:
                v_bar = 0.06067

            for node in self.solver.fluid_solver.main_model_part.GetSubModelPart("Inlet2D_Inlet").Nodes:

                vel = Vector(3)
                vel[0] = 4*v_bar*node.Y*(1-node.Y)
                vel[1] = 0.0
                vel[2] = 0.0

                node.SetSolutionStepValue(VELOCITY,0,vel)
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)

            self.solver.SetTimeStep(step)

            self.structure_main_model_part.CloneTimeStep(time)
            self.fluid_main_model_part.CloneTimeStep(time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()

            if (self.output_post == True):
                self.gid_output_structure.ExecuteInitializeSolutionStep()
                self.gid_output_fluid.ExecuteInitializeSolutionStep()

            (self.solver).Solve()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()

            if (self.output_post == True):
                self.gid_output_structure.ExecuteFinalizeSolutionStep()
                self.gid_output_fluid.ExecuteFinalizeSolutionStep()

            #TODO: decide if it shall be done only when output is processed or not
            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            if (self.output_post == True):
                if self.gid_output_structure.IsOutputStep():
                    self.gid_output_structure.PrintOutput()
                if self.gid_output_fluid.IsOutputStep():
                    self.gid_output_fluid.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

            out = out + Dt

        for process in self.list_of_processes:
            process.ExecuteFinalize()

        if (self.output_post == True):
            self.gid_output_structure.ExecuteFinalize()
            self.gid_output_fluid.ExecuteFinalize()
