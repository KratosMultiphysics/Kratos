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

# Additional imports
import timer_factory as timer_factory
import optimization_data_logger_factory as optimization_data_logger_factory

# ==============================================================================
class AlgorithmPenalizedProjection( OptimizationAlgorithm ) :

    # --------------------------------------------------------------------------
    def __init__( self, DesignSurface, DampingRegions, Analyzer, Mapper, Communicator, DataLogger, OptimizationSettings ):

        self.DesignSurface = DesignSurface
        self.Analyzer = Analyzer
        self.Mapper = Mapper
        self.Communicator = Communicator
        self.DataLogger = DataLogger        
        self.OptimizationSettings = OptimizationSettings

        self.onlyObjective = OptimizationSettings["objectives"][0]["identifier"].GetString()
        self.onlyConstraint = OptimizationSettings["constraints"][0]["identifier"].GetString()
        self.typeOfOnlyConstraint = OptimizationSettings["constraints"][0]["type"].GetString()
        self.maxIterations = OptimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1
        self.initialCorrectionScaling = OptimizationSettings["optimization_algorithm"]["correction_scaling"].GetDouble()
        self.initialStepSize = OptimizationSettings["line_search"]["step_size"].GetDouble()
        self.performDamping = OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()

        self.geometryTools = GeometryUtilities( DesignSurface )
        self.optimizationTools = OptimizationUtilities( DesignSurface, OptimizationSettings )
        if self.performDamping:        
            self.dampingUtilities = DampingUtilities( DesignSurface, DampingRegions, self.OptimizationSettings )        

    # --------------------------------------------------------------------------
    def execute( self ):
        self.__initializeOptimizationLoop()
        self.__runOptimizationLoop()
        self.__finalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeOptimizationLoop( self ):
        self.DataLogger.StartTimer()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def __runOptimizationLoop( self ):

        for optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",self.DataLogger.GetTimeStamp(),": Starting optimization iteration ", optimizationIteration)
            print(">===================================================================\n")

            self.__callCoumminicatorToCreateNewRequests()

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

            self.__updateShape()

    # --------------------------------------------------------------------------
    def __finalizeOptimizationLoop( self ):
        self.DataLogger.FinalizeDataLogging()

    # --------------------------------------------------------------------------
    def __callCoumminicatorToCreateNewRequests( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestFunctionValueOf( self.onlyObjective )
        self.Communicator.requestFunctionValueOf( self.onlyConstraint )
        self.Communicator.requestGradientOf( self.onlyObjective )
        self.Communicator.requestGradientOf( self.onlyConstraint )

    # --------------------------------------------------------------------------
    def __callAnalyzerToPerformRequestedAnalyses( self, optimizationIteration ):
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, optimizationIteration, self.Communicator )

    # --------------------------------------------------------------------------
    def __storeResultOfSensitivityAnalysisOnNodes( self ):
        gradientOfObjectiveFunction = self.Communicator.getReportedGradientOf ( self.onlyObjective )
        gradientOfConstraintFunction = self.Communicator.getReportedGradientOf ( self.onlyConstraint )
        self.__storeGradientOnNodalVariable( gradientOfObjectiveFunction, OBJECTIVE_SENSITIVITY )
        self.__storeGradientOnNodalVariable( gradientOfConstraintFunction, CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __storeGradientOnNodalVariable( self, gradients, variable_name ):
        for nodeId in gradients:
            gradient = Vector(3)
            gradient[0] = gradients[nodeId][0]
            gradient[1] = gradients[nodeId][1]
            gradient[2] = gradients[nodeId][2]
            self.DesignSurface.Nodes[nodeId].SetSolutionStepValue(variable_name,0,gradient)

    # --------------------------------------------------------------------------
    def __alignSensitivitiesToLocalSurfaceNormal( self ):
            self.geometryTools.ComputeUnitSurfaceNormals()
            self.geometryTools.ProjectNodalVariableOnUnitSurfaceNormals( OBJECTIVE_SENSITIVITY )
            self.geometryTools.ProjectNodalVariableOnUnitSurfaceNormals( CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __dampSensitivities( self ):
        self.dampingUtilities.DampNodalVariable( OBJECTIVE_SENSITIVITY )
        self.dampingUtilities.DampNodalVariable( CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __computeShapeUpdate( self ):
        self.__mapSensitivitiesToDesignSpace()

        constraintValue = self.Communicator.getReportedFunctionValueOf( self.onlyConstraint )

        if self.__isConstraintActive( constraintValue ):
            self.optimizationTools.ComputeProjectedSearchDirection()
            self.optimizationTools.CorrectProjectedSearchDirection( constraintValue )
        else:
            self.optimizationTools.ComputeSearchDirectionSteepestDescent()

        self.optimizationTools.ComputeControlPointUpdate()
        self.__mapDesignUpdateToGeometrySpace()

    # --------------------------------------------------------------------------
    def __mapSensitivitiesToDesignSpace( self ):
        self.Mapper.MapToDesignSpace( OBJECTIVE_SENSITIVITY, MAPPED_OBJECTIVE_SENSITIVITY )
        self.Mapper.MapToDesignSpace( CONSTRAINT_SENSITIVITY, MAPPED_CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __isConstraintActive( self, constraintValue ):
        if self.typeOfOnlyConstraint == "equality":
            return True
        elif constraintValue > 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __mapDesignUpdateToGeometrySpace( self ):
        self.Mapper.MapToGeometrySpace( CONTROL_POINT_UPDATE, SHAPE_UPDATE )

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

    # --------------------------------------------------------------------------
    def __updateShape( self ):
        self.optimizationTools.UpdateControlPointChangeByInputVariable( CONTROL_POINT_UPDATE )
        self.geometryTools.UpdateShapeChangeByInputVariable( SHAPE_UPDATE )

# ==============================================================================
