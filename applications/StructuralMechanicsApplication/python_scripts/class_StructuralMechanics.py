from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication


class ClassStructuralMechanics(object):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, project_parameter_file_name):
        # Read the ProjectParameters
        with open(project_parameter_file_name,'r') as parameter_file:
            self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

    def RunMainTemporalLoop():
        while self.time < self.end_time:
            self.InitializeTimeStep()
            self.SolveTimeStep()
            self.FinalizeTimeStep()

    #### Public functions defining the Interface to the CoSimulationApplication ####
    def Initialize(self):
        self.ImportAndCreateSolver()
        self.InitializeIO()
        self.ExecuteInitialize()
        self.ExecuteBeforeSolutionLoop()

    def InitializeTimeStep(self):
        self.ExecuteInitializeSolutionStep()

    def SolveTimeStep(self):
        self.SolveSolutionStep();

    def FinalizeTimeStep(self):
        self.ExecuteFinalizeSolutionStep()

    def Finalize(self):
        self.ExecuteFinalize()


    ###########################################################################
    def ImportAndCreateSolver(self, external_model_part=None):
        if external_model_part != None:
            # This is a temporary solution until the importing of the ModelPart
            # is removed from the solver (needed e.g. for Optimization)
            if (type(external_model_part) != KratosMultiphysics.ModelPart):
                raise Exception("Input is expected to be provided as a Kratos ModelPart object")
            using_external_model_part = True
        else:
            using_external_model_part = False

        ## Get echo level and parallel type
        self.echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

        ## Import parallel modules if needed
        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.mpi as KratosMPI
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication

        ## Structure model part definition
        if using_external_model_part:
            self.main_model_part = external_model_part
        else:
            main_model_part_name = ProjectParameters["problem_data"]["model_part_name"].GetString()
            self.main_model_part = ModelPart(main_model_part_name)
            self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

        ## Solver construction
        import python_solvers_wrapper_structural
        self.solver = python_solvers_wrapper_structural.CreateSolver(main_model_part, ProjectParameters)

        if not using_external_model_part:
            self.solver.AddVariables()

            ## Read the model - note that SetBufferSize is done here
            self.solver.ImportModelPart()

            ## Add AddDofs
            self.solver.AddDofs()

        ## Creation of the Kratos model (build sub_model_parts or submeshes)
        self.structure_model = KratosMultiphysics.Model()
        self.structure_model.AddModelPart(self.main_model_part)

    def InitializeIO(self):
        ## Initialize GiD  I/O
        self.output_post  = self.ProjectParameters.Has("output_configuration")
        if (self.output_post == True):
            if (self.parallel_type == "OpenMP"):
                from gid_output_process import GiDOutputProcess as output_process
            elif (self.parallel_type == "MPI"):
                from gid_output_process_mpi import GiDOutputProcessMPI as output_process

            self.gid_output = output_process(self.solver.GetComputingModelPart(),
                                             self.ProjectParameters["problem_data"]["problem_name"].GetString(),
                                             self.ProjectParameters["output_configuration"])

            self.gid_output.ExecuteInitialize()

    def ExecuteInitialize(self):
        ## Print model_part and properties
        if ((self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0)) and (self.echo_level > 1):
            KratosMultiphysics.Logger.PrintInfo("ModelPart", main_model_part)
            for properties in self.main_model_part.Properties:
                KratosMultiphysics.Logger.PrintInfo("Property " + str(properties.Id), properties)

        ## Processes construction
        import process_factory
        self.list_of_processes = process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["constraints_process_list"])
        self.list_of_processes += process_factory.KratosProcessFactory(structure_model).ConstructListOfProcesses(self.ProjectParameters["loads_process_list"])
        if (self.ProjectParameters.Has("list_other_processes") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])
        if (self.ProjectParameters.Has("json_output_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["json_output_process"])
        # Processes for tests
        if (self.ProjectParameters.Has("json_check_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["json_check_process"])
        if (self.ProjectParameters.Has("check_analytic_results_process") == True):
            self.list_of_processes += process_factory.KratosProcessFactory(self.structure_model).ConstructListOfProcesses(self.ProjectParameters["check_analytic_results_process"])

        if ((self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0)) and (self.echo_level > 1):
            count = 0
            for process in self.list_of_processes:
                count += 1
                KratosMultiphysics.Logger.PrintInfo("Process " + str(count), process)

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()

        ## Solver initialization
        self.solver.Initialize()
        self.solver.SetEchoLevel(echo_level)

    def ExecuteBeforeSolutionLoop():
        if (self.output_post == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Writing the full ProjectParameters file before solving
        if ((self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0)) and (self.echo_level > 0):
            f = open("ProjectParametersOutput.json", 'w')
            f.write(self.ProjectParameters.PrettyPrintJsonString())
            f.close()

        ## Stepping and time settings
        self.delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
        start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        self.end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        if self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True:
            self.time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = start_time
            self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0

        if (parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
            KratosMultiphysics.Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -START- ")

    def ExecuteInitializeSolutionStep(self):
        self.time += self.delta_time
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(self.time)

        if (self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
            KratosMultiphysics.Logger.PrintInfo("STEP: ", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
            KratosMultiphysics.Logger.PrintInfo("TIME: ", self.time)

        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        if (self.output_post == True):
            self.gid_output.ExecuteInitializeSolutionStep()

    def ExecuteBeforeSolve(self):
        pass

    def SolveSolutionStep(self):
        self.ExecuteBeforeSolve()
        self.solver.Solve()
        self.ExecuteAfterSolve()

    def ExecuteAfterSolve(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalizeSolutionStep()

        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

        if (self.output_post == True) and (self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()

        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()

        self.solver.SaveRestart()

    def ExecuteFinalize(self):
        for process in self.list_of_processes:
            process.ExecuteFinalize()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalize()

        if (self.parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
            KratosMultiphysics.Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -END- ")



"""
from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication  import *
from KratosMultiphysics.StructuralMechanicsApplication  import *

## Import define_output
with open("ProjectParameters.json",'r') as parameter_file:
    ProjectParameters = Parameters(parameter_file.read())

# ## Get echo level and parallel type
# echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
# parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

# ## Import parallel modules if needed
# if (parallel_type == "MPI"):
#     from KratosMultiphysics.mpi import *
#     from KratosMultiphysics.MetisApplication import *
#     from KratosMultiphysics.TrilinosApplication import *

# ## Structure model part definition
# main_model_part_name = ProjectParameters["problem_data"]["model_part_name"].GetString()
# main_model_part = ModelPart(main_model_part_name)
# main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

# ## Solver construction
# import python_solvers_wrapper_structural
# solver = python_solvers_wrapper_structural.CreateSolver(main_model_part, ProjectParameters)

# solver.AddVariables()

# ## Read the model - note that SetBufferSize is done here
# solver.ImportModelPart()

# ## Add AddDofs
# solver.AddDofs()

# ## Initialize GiD  I/O
# output_post  = ProjectParameters.Has("output_configuration")
# if (output_post == True):
#     if (parallel_type == "OpenMP"):
#         from gid_output_process import GiDOutputProcess
#         gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
#                                       ProjectParameters["problem_data"]["problem_name"].GetString() ,
#                                       ProjectParameters["output_configuration"])
#     elif (parallel_type == "MPI"):
#         from gid_output_process_mpi import GiDOutputProcessMPI
#         gid_output = GiDOutputProcessMPI(solver.GetComputingModelPart(),
#                                          ProjectParameters["problem_data"]["problem_name"].GetString() ,
#                                          ProjectParameters["output_configuration"])

#     gid_output.ExecuteInitialize()

# ## Creation of the Kratos model (build sub_model_parts or submeshes)
# StructureModel = Model()
# StructureModel.AddModelPart(main_model_part)

# ## Print model_part and properties
# if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
#     Logger.PrintInfo("ModelPart", main_model_part)
#     for properties in main_model_part.Properties:
#         Logger.PrintInfo("Property " + str(properties.Id), properties)

# ## Processes construction
# import process_factory
# list_of_processes = process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["constraints_process_list"])
# list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["loads_process_list"])
# if (ProjectParameters.Has("list_other_processes") == True):
#     list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["list_other_processes"])
# if (ProjectParameters.Has("json_output_process") == True):
#     list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["json_output_process"])

# if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
#     count = 0
#     for process in list_of_processes:
#         count += 1
#         Logger.PrintInfo("Process " + str(count), process)

# ## Processes initialization
# for process in list_of_processes:
#     process.ExecuteInitialize()

# ## Solver initialization
# solver.Initialize()
# solver.SetEchoLevel(echo_level)

# if (output_post == True):
#     gid_output.ExecuteBeforeSolutionLoop()

# for process in list_of_processes:
#     process.ExecuteBeforeSolutionLoop()

# ## Writing the full ProjectParameters file before solving
# if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 0):
#     f = open("ProjectParametersOutput.json", 'w')
#     f.write(ProjectParameters.PrettyPrintJsonString())
#     f.close()

# ## Stepping and time settings
# delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
# start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
# end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

# time = start_time
# main_model_part.ProcessInfo[STEP] = 0

# if (parallel_type == "OpenMP") or (mpi.rank == 0):
#     Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -START- ")

# Solving the problem (time integration)
while(time <= end_time):

    # time = time + delta_time
    # main_model_part.ProcessInfo[STEP] += 1
    # main_model_part.CloneTimeStep(time)

    # if (parallel_type == "OpenMP") or (mpi.rank == 0):
    #     Logger.PrintInfo("STEP: ", main_model_part.ProcessInfo[STEP])
    #     Logger.PrintInfo("TIME: ", time)

    # for process in list_of_processes:
    #     process.ExecuteInitializeSolutionStep()

    # if (output_post == True):
    #     gid_output.ExecuteInitializeSolutionStep()

    # solver.Solve()

    # for process in list_of_processes:
    #     process.ExecuteFinalizeSolutionStep()

    # if (output_post == True):
    #     gid_output.ExecuteFinalizeSolutionStep()

    # for process in list_of_processes:
    #     process.ExecuteBeforeOutputStep()

    # if (output_post == True) and (gid_output.IsOutputStep()):
    #     gid_output.PrintOutput()

    # for process in list_of_processes:
    #     process.ExecuteAfterOutputStep()

# for process in list_of_processes:
#     process.ExecuteFinalize()

# if (output_post == True):
#     gid_output.ExecuteFinalize()

# if (parallel_type == "OpenMP") or (mpi.rank == 0):
#     Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -END- ")


"""