from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

#import define_output
parameter_file = open("ProjectParameters_with_mass_constraint.json",'r')
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

# Create solver for handling mesh-motion
mesh_solver_module = __import__(ProjectParameters["mesh_solver_settings"]["solver_type"].GetString())
mesh_solver = mesh_solver_module.CreateSolver(main_model_part, ProjectParameters["mesh_solver_settings"])
mesh_solver.AddVariables()

# Create solver to perform structural analysis
solver_module = __import__(ProjectParameters["structure_solver_settings"]["solver_type"].GetString())
CSM_solver = solver_module.CreateSolver(main_model_part, ProjectParameters["structure_solver_settings"])
CSM_solver.AddVariables()
CSM_solver.ImportModelPart()

# Add degrees of freedom
CSM_solver.AddDofs()
mesh_solver.AddDofs()

# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["structure_solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["structure_solver_settings"]["processes_sub_model_part_list"][i].GetString()
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

        mesh_solver.Initialize()
        mesh_solver.SetEchoLevel(echo_level)        

        for responseFunctionId in listOfResponseFunctions:
            listOfResponseFunctions[responseFunctionId].Initialize()

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

            print("\n> Starting ALEApplication to update the mesh")
            startTime = timer.time()
            self.updateMeshForAnalysis()
            print("> Time needed for updating the mesh = ",round(timer.time() - startTime,2),"s")

            print("\n> Starting StructuralMechanicsApplication to solve structure")
            startTime = timer.time()
            self.solveStructure( optimizationIteration )
            print("> Time needed for solving the structure = ",round(timer.time() - startTime,2),"s")

            print("\n> Starting calculation of strain energy")
            startTime = timer.time()                    
            listOfResponseFunctions["strain_energy"].CalculateValue()
            print("> Time needed for calculation of strain energy = ",round(timer.time() - startTime,2),"s")

            communicator.reportFunctionValue("strain_energy", listOfResponseFunctions["strain_energy"].GetValue())  

        # Calculation of value of constraint function
        if communicator.isRequestingFunctionValueOf("mass"):

            print("\n> Starting calculation of value of mass constraint")
            listOfResponseFunctions["mass"].CalculateValue()
            constraintFunctionValue = listOfResponseFunctions["mass"].GetValue() - listOfResponseFunctions["mass"].GetInitialValue()
            print("> Time needed for calculation of value of mass constraint = ",round(timer.time() - startTime,2),"s")

            communicator.reportFunctionValue("mass", constraintFunctionValue)
            communicator.setFunctionReferenceValue("mass", listOfResponseFunctions["mass"].GetInitialValue())

        # Calculation of gradients of objective function
        if communicator.isRequestingGradientOf("strain_energy"): 

            print("\n> Starting calculation of gradient of objective function")
            startTime = timer.time()               
            listOfResponseFunctions["strain_energy"].CalculateGradient()
            print("> Time needed for calculating gradient of objective function = ",round(timer.time() - startTime,2),"s")
            
            gradientForCompleteModelPart = listOfResponseFunctions["strain_energy"].GetGradient()
            gradientOnDesignSurface = {}
            for node in currentDesign.Nodes:
                gradientOnDesignSurface[node.Id] = gradientForCompleteModelPart[node.Id]

            # If contribution from mesh-motion to gradient shall be considered
            # self.computeAndAddMeshDerivativesToGradient(gradientOnDesignSurface, gradientForCompleteModelPart)                 

            communicator.reportGradient("strain_energy", gradientOnDesignSurface)   

        # Calculation of gradients of constraint function  
        if communicator.isRequestingGradientOf("mass"):       

            print("\n> Starting calculation of gradient of constraint function")
            startTime = timer.time()    
            listOfResponseFunctions["mass"].CalculateGradient()
            print("> Time needed for calculating gradient of constraint function = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = listOfResponseFunctions["mass"].GetGradient()
            gradientOnDesignSurface = {}
            for node in currentDesign.Nodes:
                gradientOnDesignSurface[node.Id] = gradientForCompleteModelPart[node.Id]

            communicator.reportGradient("mass", gradientOnDesignSurface) 

    # --------------------------------------------------------------------------
    def initializeNewSolutionStep( self, optimizationIteration ):
        main_model_part.CloneTimeStep( optimizationIteration )

    # --------------------------------------------------------------------------
    def updateMeshForAnalysis( self ):

        # Extract surface nodes
        sub_model_part_name = "surface_nodes"     
        GeometryUtilities(main_model_part).ExtractSurfaceNodes(sub_model_part_name)

        # Apply shape update as boundary condition for computation of mesh displacement 
        for node in main_model_part.GetSubModelPart(sub_model_part_name).Nodes:
            node.Fix(MESH_DISPLACEMENT_X)
            node.Fix(MESH_DISPLACEMENT_Y)
            node.Fix(MESH_DISPLACEMENT_Z)              
            disp = Vector(3)
            disp[0] = node.GetSolutionStepValue(SHAPE_UPDATE_X)
            disp[1] = node.GetSolutionStepValue(SHAPE_UPDATE_Y)
            disp[2] = node.GetSolutionStepValue(SHAPE_UPDATE_Z)
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,disp)

        # Solve for mesh-update
        mesh_solver.Solve()

        # Update reference mesh (Since shape updates are imposed as incremental quantities)
        mesh_solver.get_mesh_motion_solver().UpdateReferenceMesh()

        # Log absolute mesh displacement
        for node in main_model_part.Nodes:
            mesh_change = Vector(3)
            mesh_change[0] = node.GetSolutionStepValue(MESH_CHANGE_X) + node.GetSolutionStepValue(MESH_DISPLACEMENT_X)
            mesh_change[1] = node.GetSolutionStepValue(MESH_CHANGE_Y) + node.GetSolutionStepValue(MESH_DISPLACEMENT_Y)
            mesh_change[2] = node.GetSolutionStepValue(MESH_CHANGE_Z) + node.GetSolutionStepValue(MESH_DISPLACEMENT_Z)
            node.SetSolutionStepValue(MESH_CHANGE,0,mesh_change)                                      

    # --------------------------------------------------------------------------
    def solveStructure( self, optimizationIteration ): 

        # processes to be executed at the begining of the solution step
        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

        self.gid_output.ExecuteInitializeSolutionStep()
            
        # Actual solution
        CSM_solver.Solve()

        # processes to be executed at the end of the solution step
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()
        
        # write output results GiD: (frequency writing is controlled internally)
        if(self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()

        self.gid_output.ExecuteFinalizeSolutionStep()
                        
        # processes to be executed after witting the output
        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()            

    # --------------------------------------------------------------------------
    def computeAndAddMeshDerivativesToGradient( self, gradientOnDesignSurface, gradientForCompleteModelPart):
    
        # Here we solve the pseudo-elastic mesh-motion system again using modified BCs
        # The contributions from the mesh derivatives appear as reaction forces
        for node in main_model_part.Nodes:

            # Apply dirichlet conditions
            if node.Id in gradientOnDesignSurface.keys():
                node.Fix(MESH_DISPLACEMENT_X)
                node.Fix(MESH_DISPLACEMENT_Y)
                node.Fix(MESH_DISPLACEMENT_Z)              
                xs = Vector(3)
                xs[0] = 0.0
                xs[1] = 0.0
                xs[2] = 0.0
                node.SetSolutionStepValue(MESH_DISPLACEMENT,0,xs)
            # Apply RHS conditions       
            else:
                rhs = Vector(3)
                rhs[0] = gradientForCompleteModelPart[node.Id][0]
                rhs[1] = gradientForCompleteModelPart[node.Id][1]
                rhs[2] = gradientForCompleteModelPart[node.Id][2]
                node.SetSolutionStepValue(MESH_RHS,0,rhs)

        # Solve mesh-motion problem with previously modified BCs
        mesh_solver.Solve()

        # Compute and add gradient contribution from mesh motion
        for node_id in gradientOnDesignSurface.keys():
            node = main_model_part.Nodes[node_id]
            sens_contribution = Vector(3)
            sens_contribution = node.GetSolutionStepValue(MESH_REACTION)
            gradientOnDesignSurface[node.Id] = gradientOnDesignSurface[node_id] + sens_contribution    
    
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