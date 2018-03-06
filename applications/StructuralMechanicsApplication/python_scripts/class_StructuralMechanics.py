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
        ## Import define_output
        with open(project_parameter_file_name,'r') as parameter_file:
            ProjectParameters = Parameters(parameter_file.read())

    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

    def RunMainTemporalLoop():
        while time < end_time:
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
        echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
        parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

        ## Import parallel modules if needed
        if (parallel_type == "MPI"):
            from KratosMultiphysics.mpi import *
            from KratosMultiphysics.MetisApplication import *
            from KratosMultiphysics.TrilinosApplication import *

        ## Structure model part definition
        if using_external_model_part:
            main_model_part = external_model_part
        else:
            main_model_part_name = ProjectParameters["problem_data"]["model_part_name"].GetString()
            main_model_part = ModelPart(main_model_part_name)
            main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

        ## Solver construction
        import python_solvers_wrapper_structural
        solver = python_solvers_wrapper_structural.CreateSolver(main_model_part, ProjectParameters)

        if not using_external_model_part:
            solver.AddVariables()

            ## Read the model - note that SetBufferSize is done here
            solver.ImportModelPart()

            ## Add AddDofs
            solver.AddDofs()

        ## Creation of the Kratos model (build sub_model_parts or submeshes)
        StructureModel = Model()
        StructureModel.AddModelPart(main_model_part)

    def InitializeIO(self):
        ## Initialize GiD  I/O
        output_post  = ProjectParameters.Has("output_configuration")
        if (output_post == True):
            if (parallel_type == "OpenMP"):
                from gid_output_process import GiDOutputProcess
                gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                                            ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                            ProjectParameters["output_configuration"])
            elif (parallel_type == "MPI"):
                from gid_output_process_mpi import GiDOutputProcessMPI
                gid_output = GiDOutputProcessMPI(solver.GetComputingModelPart(),
                                                ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                                ProjectParameters["output_configuration"])

            gid_output.ExecuteInitialize()

    def ExecuteInitialize(self):
        ## Print model_part and properties
        if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
            Logger.PrintInfo("ModelPart", main_model_part)
            for properties in main_model_part.Properties:
                Logger.PrintInfo("Property " + str(properties.Id), properties)

        ## Processes construction
        import process_factory
        list_of_processes = process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["constraints_process_list"])
        list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["loads_process_list"])
        if (ProjectParameters.Has("list_other_processes") == True):
            list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["list_other_processes"])
        if (ProjectParameters.Has("json_output_process") == True):
            list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["json_output_process"])

        if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
            count = 0
            for process in list_of_processes:
                count += 1
                Logger.PrintInfo("Process " + str(count), process)

        ## Processes initialization
        for process in list_of_processes:
            process.ExecuteInitialize()

        ## Solver initialization
        solver.Initialize()
        solver.SetEchoLevel(echo_level)

    def ExecuteBeforeSolutionLoop():
        if (output_post == True):
            gid_output.ExecuteBeforeSolutionLoop()

        for process in list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Writing the full ProjectParameters file before solving
        if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 0):
            f = open("ProjectParametersOutput.json", 'w')
            f.write(ProjectParameters.PrettyPrintJsonString())
            f.close()

        ## Stepping and time settings
        delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
        start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
        end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = start_time
        main_model_part.ProcessInfo[STEP] = 0

        if (parallel_type == "OpenMP") or (mpi.rank == 0):
            Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -START- ")

    def ExecuteInitializeSolutionStep(self):
        time = time + delta_time
        main_model_part.ProcessInfo[STEP] += 1
        main_model_part.CloneTimeStep(time)

        if (parallel_type == "OpenMP") or (mpi.rank == 0):
            Logger.PrintInfo("STEP: ", main_model_part.ProcessInfo[STEP])
            Logger.PrintInfo("TIME: ", time)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        if (output_post == True):
            gid_output.ExecuteInitializeSolutionStep()

    def SolveSolutionStep(self):
        self.ExecuteBeforeSolve()
        self.solver.Solve()
        self.ExecuteAfterSolve()

    def ExecuteFinalizeSolutionStep(self):
        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        if (output_post == True):
            gid_output.ExecuteFinalizeSolutionStep()

        for process in list_of_processes:
            process.ExecuteBeforeOutputStep()

        if (output_post == True) and (gid_output.IsOutputStep()):
            gid_output.PrintOutput()

        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        for process in list_of_processes:
            process.ExecuteFinalize()

        if (output_post == True):
            gid_output.ExecuteFinalize()

        if (parallel_type == "OpenMP") or (mpi.rank == 0):
            Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -END- ")



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