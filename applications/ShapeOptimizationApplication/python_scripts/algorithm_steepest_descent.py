# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Import algorithm base classes
from algorithm_base import OptimizationAlgorithm

# ==============================================================================
class AlgorithmSteepestDescent( OptimizationAlgorithm ) :

    # --------------------------------------------------------------------------
    def __init__( self, 
                  DesignSurface, 
                  DampingRegions, 
                  Analyzer, 
                  Communicator, 
                  MeshController, 
                  Mapper, 
                  DataLogger, 
                  OptimizationSettings ):

        self.DesignSurface = DesignSurface
        self.Analyzer = Analyzer
        self.Communicator = Communicator
        self.MeshController = MeshController
        self.Mapper = Mapper
        self.DataLogger = DataLogger
        self.OptimizationSettings = OptimizationSettings

        self.maxIterations = OptimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1
        self.onlyObjective = OptimizationSettings["objectives"][0]["identifier"].GetString()
        self.initialStepSize = OptimizationSettings["line_search"]["step_size"].GetDouble()
        self.performDamping = OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()

        self.GeometryUtilities = GeometryUtilities( DesignSurface )
        self.OptimizationUtilities = OptimizationUtilities( DesignSurface, OptimizationSettings )
        if self.performDamping:
            self.dampingUtilities = DampingUtilities( DesignSurface, DampingRegions, self.OptimizationSettings )
            
    # --------------------------------------------------------------------------
    def execute( self ):
        self.__initializeOptimizationLoop()
        self.__runOptimizationLoop()
        self.__finalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeOptimizationLoop( self ):
        self.MeshController.Initialize()
        self.DataLogger.StartTimer()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def __runOptimizationLoop( self ):

        for optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",self.DataLogger.GetTimeStamp(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            self.__updateMeshAccordingCurrentShapeUpdate()    

            self.__callCommunicatorToRequestNewAnalyses()

            self.__callAnalyzerToPerformRequestedAnalyses( optimizationIteration )

            self.__storeResultOfSensitivityAnalysisOnNodes()

            self.__alignSensitivitiesToLocalSurfaceNormal()

            if self.performDamping:
                self.__dampSensitivities()

            self.__computeShapeUpdate()

            if self.performDamping:
                self.__dampShapeUpdate()

            self.__logCurrentOptimizationStep( optimizationIteration )

            self.__timeOptimizationStep()

            if self.__isAlgorithmConverged( optimizationIteration ):
                break

    # --------------------------------------------------------------------------
    def __finalizeOptimizationLoop( self ):
        self.DataLogger.FinalizeDataLogging()
        self.MeshController.Finalize()

    # --------------------------------------------------------------------------
    def __updateMeshAccordingCurrentShapeUpdate( self ):
        self.MeshController.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )  

    # --------------------------------------------------------------------------
    def __callCommunicatorToRequestNewAnalyses( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestFunctionValueOf( self.onlyObjective )
        self.Communicator.requestGradientOf( self.onlyObjective )

    # --------------------------------------------------------------------------
    def __callAnalyzerToPerformRequestedAnalyses( self, optimizationIteration ):
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, optimizationIteration, self.Communicator )
        self.__ResetPossibleMeshDisplacementDuringAnalysis()

    # --------------------------------------------------------------------------
    def __ResetPossibleMeshDisplacementDuringAnalysis( self ):
        self.MeshController.ResetMeshDisplacement()

    # --------------------------------------------------------------------------
    def __storeResultOfSensitivityAnalysisOnNodes( self ):
        gradientOfObjectiveFunction = self.Communicator.getReportedGradientOf ( self.onlyObjective )
        self.__storeGradientOnNodalVariable( gradientOfObjectiveFunction, OBJECTIVE_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __storeGradientOnNodalVariable( self , givenGradient, variable_name ):
        for nodeId in givenGradient:
            gradient = Vector(3)
            gradient[0] = givenGradient[nodeId][0]
            gradient[1] = givenGradient[nodeId][1]
            gradient[2] = givenGradient[nodeId][2]
            self.DesignSurface.Nodes[nodeId].SetSolutionStepValue(variable_name,0,gradient)

    # --------------------------------------------------------------------------
    def __alignSensitivitiesToLocalSurfaceNormal( self ):
            self.GeometryUtilities.ComputeUnitSurfaceNormals()
            self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( OBJECTIVE_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __dampSensitivities( self ):
        self.dampingUtilities.DampNodalVariable( OBJECTIVE_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __computeShapeUpdate( self ):
        self.__mapSensitivitiesToDesignSpace()
        self.OptimizationUtilities.ComputeSearchDirectionSteepestDescent()
        self.OptimizationUtilities.ComputeControlPointUpdate()
        self.__mapDesignUpdateToGeometrySpace()
        self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def __mapSensitivitiesToDesignSpace( self ):
        self.Mapper.MapToDesignSpace( OBJECTIVE_SENSITIVITY, MAPPED_OBJECTIVE_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __mapDesignUpdateToGeometrySpace( self ):
        self.Mapper.MapToGeometrySpace( CONTROL_POINT_UPDATE, SHAPE_UPDATE )

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges( self ):
        self.OptimizationUtilities.AddFirstVariableToSecondVariable( CONTROL_POINT_UPDATE, CONTROL_POINT_CHANGE )        
        self.OptimizationUtilities.AddFirstVariableToSecondVariable( SHAPE_UPDATE, SHAPE_CHANGE )

    # --------------------------------------------------------------------------
    def __dampShapeUpdate( self ):
        self.dampingUtilities.DampNodalVariable( SHAPE_UPDATE )

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep( self, optimizationIteration ):
        self.DataLogger.LogCurrentData( optimizationIteration )

    # --------------------------------------------------------------------------
    def __timeOptimizationStep( self ):
        print("\n> Time needed for current optimization step = ", self.DataLogger.GetLapTime(), "s")
        print("> Time needed for total optimization so far = ", self.DataLogger.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged( self, optimizationIteration ):

        if optimizationIteration > 1 :

            # Check if maximum iterations were reached
            if optimizationIteration == self.maxIterations:
                print("\n> Maximal iterations of optimization problem reached!")
                return True

            relativeChangeOfObjectiveValue = self.DataLogger.GetValue( "RELATIVE_CHANGE_OF_OBJECTIVE_VALUE" )

            # Check for relative tolerance
            relativeTolerance = self.OptimizationSettings["optimization_algorithm"]["relative_tolerance"].GetDouble()
            if abs(relativeChangeOfObjectiveValue) < relativeTolerance:
                print("\n> Optimization problem converged within a relative objective tolerance of ",relativeTolerance,"%.")
                return True

            # Check if value of objective increases
            if relativeChangeOfObjectiveValue > 0:
                print("\n> Value of objective function increased!")
                return False          

# ==============================================================================
