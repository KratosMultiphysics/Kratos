from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.ChimeraApplication as KratosChimera
from KratosMultiphysics.ExternalSolversApplication import *

import os
import KratosMultiphysics.KratosUnittest as KratosUnittest

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TestChimeraMonolithicSimple(KratosUnittest.TestCase):

    def test_MonolithicSimple(self):
        self._setUp()
        self._runTest()
        self._tearDown()
        self._checkResults()

    def _setUp(self):
        self.check_tolerance = 1e-6
        self.work_folder = "chimera_monolithic_simple_test"
        self.reference_file = "reference_chimera_monolithic_simple_test"
        self.print_reference_values = False


    def _runTest(self):
        with WorkFolderScope(self.work_folder):
            parameter_file = open("test_chimera_monolithic_simple_ProjectParameters.json",'r')
            self.ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())
            ProjectParameters = self.ProjectParameters

            ## Get echo level and parallel type
            echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
            parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()
            # Import KratosMPI if needed

            solver_settings = ProjectParameters["solver_settings"]
            if not solver_settings.Has("domain_size"):
                KratosMultiphysics.Logger.PrintInfo("FluidDynamicsAnalysis", "Using the old way to pass the domain_size, this will be removed!")
                solver_settings.AddEmptyValue("domain_size")
                solver_settings["domain_size"].SetInt(ProjectParameters["problem_data"]["domain_size"].GetInt())

            if not solver_settings.Has("model_part_name"):
                KratosMultiphysics.Logger.PrintInfo("FluidDynamicsAnalysis", "Using the old way to pass the model_part_name, this will be removed!")
                solver_settings.AddEmptyValue("model_part_name")
                solver_settings["model_part_name"].SetString(ProjectParameters["problem_data"]["model_part_name"].GetString())

            fluid_model = KratosMultiphysics.Model()


            ## Solver construction
            import python_solvers_wrapper_fluid_chimera
            solver = python_solvers_wrapper_fluid_chimera.CreateSolver(fluid_model, ProjectParameters)

            solver.AddVariables()

            ## Read the model - note that SetBufferSize is done here
            solver.ImportModelPart()
            solver.PrepareModelPart()

            main_model_part = fluid_model[ProjectParameters["problem_data"]["model_part_name"].GetString()]
            self.main_model_part = main_model_part
            ## Add AddDofs
            solver.AddDofs()

            ## Modelparts
            background = main_model_part.GetSubModelPart("GENERIC_background")
            patch= main_model_part.GetSubModelPart("GENERIC_patch")
            patchBoundary = main_model_part.GetSubModelPart("GENERIC_patchBoundary")
            structure = main_model_part.GetSubModelPart("NoSlip2D_structure")

            ## Processes construction
            import process_factory
            # "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
            # Note 1: gravity is firstly constructed. Outlet process might need its information.
            # Note 2: conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
            list_of_processes =  process_factory.KratosProcessFactory(fluid_model).ConstructListOfProcesses( ProjectParameters["gravity"] )
            list_of_processes += process_factory.KratosProcessFactory(fluid_model).ConstructListOfProcesses( ProjectParameters["initial_conditions_process_list"] )
            list_of_processes += process_factory.KratosProcessFactory(fluid_model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )
            list_of_processes += process_factory.KratosProcessFactory(fluid_model).ConstructListOfProcesses( ProjectParameters["auxiliar_process_list"] )

            if (echo_level > 1) and ((parallel_type == "OpenMP") or (mpi.rank == 0)):
                for process in list_of_processes:
                    print(process)

            ## Processes initialization
            for process in list_of_processes:
                process.ExecuteInitialize()

            ## Solver initialization
            solver.Initialize()

            ## Stepping and time settings
            start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
            #Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
            end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

            time = start_time
            step = 0

            for process in list_of_processes:
                process.ExecuteBeforeSolutionLoop()

            ChimeraParamaters = KratosMultiphysics.Parameters(R"""{
                                "process_name":"chimera",
                                "Chimera_levels" : [ [{ "model_part_name":"GENERIC_background",
                                                        "model_part_boundary_name":"GENERIC_DomainBoundary",
                                                        "model_part_inside_boundary_name" :"GENERIC_DomainBoundary"
		                        }],
                                       [{ "model_part_name":"GENERIC_patch",
                                        "model_part_boundary_name":"GENERIC_patchBoundary",
                                        "model_part_inside_boundary_name":"GENERIC_StrctureOne"
                                        }]
                                        ],
                                "type" : "nearest_element",
                                "IsWeak" : true,
                                "pressure_coupling" : "all",
                                "pressure_coupling_node" : 0.0,
                                "overlap_distance":0.01
                                }""")

            ChimeraProcess = KratosChimera.ApplyChimeraProcessMonolithic2d(main_model_part,ProjectParameters["chimera"][0])

            while(time <= end_time):

                time = solver.AdvanceInTime(time)

                ChimeraProcess.ExecuteInitializeSolutionStep()

                for process in list_of_processes:
                    process.ExecuteInitializeSolutionStep()

                solver.Solve()

                for process in list_of_processes:
                    process.ExecuteFinalizeSolutionStep()

                ChimeraProcess.ExecuteFinalizeSolutionStep()

                for process in list_of_processes:
                    process.ExecuteBeforeOutputStep()

                for process in list_of_processes:
                    process.ExecuteAfterOutputStep()

            for process in list_of_processes:
                process.ExecuteFinalize()

    def _tearDown(self):
        with WorkFolderScope(self.work_folder):
            try:
                os.remove(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')
            except FileNotFoundError as e:
                pass

    def _checkResults(self):
        with WorkFolderScope(self.work_folder):
            if self.print_reference_values:
                with open(self.reference_file+'.csv','w') as ref_file:
                    ref_file.write("#ID, VELOCITY_X, VELOCITY_Y\n")
                    for node in self.main_model_part.Nodes:
                        vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                        ref_file.write("{0}, {1}, {2}\n".format(node.Id, vel[0], vel[1]))
            else:
                with open(self.reference_file+'.csv','r') as reference_file:
                    reference_file.readline() # skip header
                    line = reference_file.readline()

                    for node in self.main_model_part.Nodes:
                        values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                        node_id = values[0]
                        reference_vel_x = values[1]
                        reference_vel_y = values[2]

                        velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                        self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                        self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)

                        line = reference_file.readline()
                    if line != '': # If we did not reach the end of the reference file
                        self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

if __name__ == '__main__':
    KratosUnittest.main()
