from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#defining the model_part
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

# Create an optimizer
# Note that internally variables related to the optimizer are added to the model part
optimizerFactory = __import__("optimizer_factory")
optimizer = optimizerFactory.CreateOptimizer( main_model_part, ProjectParameters["optimization_settings"] )

# Create solver for all response functions specified in the optimization settings
# Note that internally variables related to the individual functions are added to the model part
responseFunctionFactory = __import__("response_function_factory")
listOfResponseFunctions = responseFunctionFactory.CreateListOfResponseFunctions( main_model_part, ProjectParameters["optimization_settings"] )

# Create solver to perform structural analysis
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
CSM_solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])
CSM_solver.AddVariables()
CSM_solver.ImportModelPart()

#add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
CSM_solver.AddDofs()

# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    if( main_model_part.HasSubModelPart(part_name) ):
        Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

# ======================================================================================================================================
# Analyzer
# ======================================================================================================================================

class kratosCSMAnalyzer( (__import__("analyzer_base")).analyzerBaseClass ):

    # --------------------------------------------------------------------------
    def __init__( self ):

        self.initializeGIDOutput()
        self.initializeProcesses()
        self.initializeSolutionLoop()

    # --------------------------------------------------------------------------
    def initializeProcesses( self ):

        import process_factory
        #the process order of execution is important
        self.list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
        self.list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )
        if(ProjectParameters.Has("problem_process_list")):
            self.list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["problem_process_list"] )
        if(ProjectParameters.Has("output_process_list")):
            self.list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["output_process_list"] )

        #print list of constructed processes
        if(echo_level>1):
            for process in self.list_of_processes:
                print(process)

        for process in self.list_of_processes:
            process.ExecuteInitialize()

    # --------------------------------------------------------------------------
    def initializeGIDOutput( self ):

        computing_model_part = CSM_solver.GetComputingModelPart()
        problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

        from gid_output_process import GiDOutputProcess
        output_settings = ProjectParameters["output_configuration"]
        self.gid_output = GiDOutputProcess(computing_model_part, problem_name, output_settings)

        self.gid_output.ExecuteInitialize()

    # --------------------------------------------------------------------------
    def initializeSolutionLoop( self ):

        ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        CSM_solver.Initialize()
        CSM_solver.SetEchoLevel(echo_level)

        for responseFunctionId in listOfResponseFunctions:
            listOfResponseFunctions[responseFunctionId].initialize()

        # Start process
        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Set results when are written in a single file
        self.gid_output.ExecuteBeforeSolutionLoop()

    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        # Calculation of value of objective function
        if communicator.isRequestingFunctionValueOf("strain_energy"):

            self.initializeNewSolutionStep( optimizationIteration )

            print("\n> Starting to update the mesh")
            startTime = timer.time()
            self.updateMeshForAnalysis( currentDesign )
            print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

            print("\n> Starting StructuralMechanicsApplication to solve structure")
            startTime = timer.time()
            self.solveStructure( optimizationIteration )
            print("> Time needed for solving the structure = ",round(timer.time() - startTime,2),"s")

            print("\n> Starting calculation of response value")
            startTime = timer.time()
            listOfResponseFunctions["strain_energy"].calculate_value()
            print("> Time needed for calculation of response value = ",round(timer.time() - startTime,2),"s")

            communicator.reportFunctionValue("strain_energy", listOfResponseFunctions["strain_energy"].get_value())

        # Calculation of gradient of objective function
        if communicator.isRequestingGradientOf("strain_energy"):

            print("\n> Starting calculation of gradients")
            startTime = timer.time()
            listOfResponseFunctions["strain_energy"].calculate_gradient()
            print("> Time needed for calculating gradients = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = listOfResponseFunctions["strain_energy"].get_gradient()
            gradientOnDesignSurface = {}
            for node in currentDesign.Nodes:
                gradientOnDesignSurface[node.Id] = gradientForCompleteModelPart[node.Id]

            communicator.reportGradient("strain_energy", gradientOnDesignSurface)

    # --------------------------------------------------------------------------
    def initializeNewSolutionStep( self, optimizationIteration ):
        main_model_part.CloneTimeStep( optimizationIteration )

    # --------------------------------------------------------------------------
    def updateMeshForAnalysis( self, currentDesign ):
        for node in currentDesign.Nodes:
            node.X0 = node.X0 + node.GetSolutionStepValue(SHAPE_UPDATE_X)
            node.Y0 = node.Y0 + node.GetSolutionStepValue(SHAPE_UPDATE_Y)
            node.Z0 = node.Z0 + node.GetSolutionStepValue(SHAPE_UPDATE_Z)

    # --------------------------------------------------------------------------
    def solveStructure( self, optimizationIteration ):

        # processes to be executed at the begining of the solution step
        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        self.gid_output.ExecuteInitializeSolutionStep()

        # Actual solution
        CSM_solver.Solve()

        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        self.gid_output.ExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

        # write output results GiD: (frequency writing is controlled internally)
        if(self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()

        # processes to be executed after witting the output
        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()

    # --------------------------------------------------------------------------
    def finalizeSolutionLoop( self ):
        for process in self.list_of_processes:
            process.ExecuteFinalize()
        self.gid_output.ExecuteFinalize()

    # --------------------------------------------------------------------------

structureAnalyzer = kratosCSMAnalyzer()

# ======================================================================================================================================
# Optimization
# ======================================================================================================================================

optimizer.importAnalyzer( structureAnalyzer )
optimizer.optimize()
structureAnalyzer.finalizeSolutionLoop()

# ======================================================================================================================================
