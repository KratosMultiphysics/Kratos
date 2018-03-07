from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# TODO replace this with the future default solver from the applications
import structural_mechanics_analysis

# For time measures
import time as timer

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

parameter_file = open("optimization_parameters.json",'r')
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
#responseFunctionFactory = __import__("response_function_factory")
#listOfResponseFunctions = responseFunctionFactory.CreateListOfResponseFunctions( main_model_part, ProjectParameters["optimization_settings"] )

responseSettings = ProjectParameters["optimization_settings"]["objectives"][0]["adjoint_response_settings"]
# Create the primal solver
ProjectParametersPrimal = Parameters( open(responseSettings["primal_settings"].GetString(),'r').read() )
primal_solver = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal)

# Create the adjoint solver
ProjectParametersAdjoint = Parameters( open(responseSettings["adjoint_settings"].GetString(),'r').read() )
adjoint_solver = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)

# ======================================================================================================================================
# Analyzer
# ======================================================================================================================================

class kratosAdjointAnalyzer( (__import__("analyzer_base")).analyzerBaseClass ):

    # --------------------------------------------------------------------------
    # NOTES:
    # Currently a separate modelpart is created for:
    #   - optimizer
    #   - primal solver
    #   - adjoint solver
    # The nodal coordinates are synchronized at each analyzer call

    # --------------------------------------------------------------------------
    def initializeBeforeOptimizationLoop( self ):
        primal_solver.Initialize()
        adjoint_solver.Initialize()

    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        # update mesh of primal and adjoint analysis
        primal_solver.main_model_part.CloneTimeStep( optimizationIteration )
        for node in primal_solver.main_model_part.Nodes:
            ref_node = main_model_part.Nodes[node.Id]
            node.X0 = ref_node.X0
            node.Y0 = ref_node.Y0
            node.Z0 = ref_node.Z0
            node.X = ref_node.X
            node.Y = ref_node.Y
            node.Z = ref_node.Z

        adjoint_solver.main_model_part.CloneTimeStep( optimizationIteration )
        for node in adjoint_solver.main_model_part.Nodes:
            ref_node = main_model_part.Nodes[node.Id]
            node.X0 = ref_node.X0
            node.Y0 = ref_node.Y0
            node.Z0 = ref_node.Z0
            node.X = ref_node.X
            node.Y = ref_node.Y
            node.Z = ref_node.Z

        # Calculation of value of strain energy
        if communicator.isRequestingValueOf("strain_energy"):

            print("\n> Starting StructuralMechanicsApplication to solve structure")
            startTime = timer.time()
            primal_solver.InitializeTimeStep()
            primal_solver.SolveTimeStep()
            primal_solver.FinalizeTimeStep()
            print("> Time needed for solving the structure = ",round(timer.time() - startTime,2),"s")

            #TODO transfer primal solution

            print("\n> Starting calculation of strain energy")
            startTime = timer.time()
            print(primal_solver.main_model_part)
            print(adjoint_solver.main_model_part)
            # TODO input process was not yet executed so here only primal modelpart has values, anyway the value should be calculated on top of the primal
            value = adjoint_solver.solver.response_function.CalculateValue(primal_solver.main_model_part)
            print("> Time needed for calculation of strain energy = ",round(timer.time() - startTime,2),"s")

            print(" OBEJCTIVE VALUE =", value)
            communicator.reportValue("strain_energy", value)

        # Calculation of gradient of strain energy
        if communicator.isRequestingGradientOf("strain_energy"):

            print("\n> Starting calculation of gradients of strain energy")
            startTime = timer.time()
            adjoint_solver.InitializeTimeStep()
            adjoint_solver.SolveTimeStep()
            adjoint_solver.FinalizeTimeStep()
            print("> Time needed for calculating gradients of strain energy = ",round(timer.time() - startTime,2),"s")

            gradientOnDesignSurface = {}
            for node in adjoint_solver.main_model_part.Nodes:
                grad = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
                gradientOnDesignSurface[node.Id] = grad

            communicator.reportGradient("strain_energy", gradientOnDesignSurface)

        MeshControllerUtilities( primal_solver.main_model_part ).SetDeformationVariablesToZero()

        MeshControllerUtilities( adjoint_solver.main_model_part ).SetDeformationVariablesToZero()
        # TODO reset also ADJOINT_DISPLACEMENT and ADJOINT_ROTATION

    # --------------------------------------------------------------------------
    def finalizeAfterOptimizationLoop( self ):
        primal_solver.Finalize()
        adjoint_solver.Finalize()

structureAnalyzer = kratosAdjointAnalyzer()

# ======================================================================================================================================
# Optimization
# ======================================================================================================================================

optimizer.importAnalyzer( structureAnalyzer )
optimizer.optimize()

# ======================================================================================================================================