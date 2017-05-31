# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
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
    def __init__( self, designSurface, dampingRegions, analyzer, mapper, communicator, optimizationSettings ):

        self.designSurface = designSurface
        self.analyzer = analyzer
        self.mapper = mapper
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings

        self.onlyObjective = optimizationSettings["objectives"][0]["identifier"].GetString()
        self.onlyConstraint = optimizationSettings["constraints"][0]["identifier"].GetString()
        self.typeOfOnlyConstraint = optimizationSettings["constraints"][0]["type"].GetString()          
        self.maxIterations = optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1        
        self.initialCorrectionScaling = optimizationSettings["optimization_algorithm"]["correction_scaling"].GetDouble()
        self.initialStepSize = optimizationSettings["line_search"]["step_size"].GetDouble()
        self.performDamping = optimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()

        self.geometryTools = GeometryUtilities( designSurface )
        self.optimizationTools = OptimizationUtilities( designSurface, optimizationSettings )
        self.dampingUtilities = DampingUtilities( designSurface, dampingRegions, self.optimizationSettings )        

        self.timer = timer_factory.CreateTimer()
        self.dataLogger = optimization_data_logger_factory.CreateDataLogger( designSurface, communicator, optimizationSettings, self.timer )

    # --------------------------------------------------------------------------
    def execute( self ):
        self.__initializeOptimizationLoop()
        self.__startOptimizationLoop()
        self.__finalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeOptimizationLoop( self ):   
        self.timer.startTimer()
        self.dataLogger.initializeDataLogging() 

    # --------------------------------------------------------------------------
    def __startOptimizationLoop( self ):

        for optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",self.timer.getTimeStamp(),": Starting optimization iteration ", optimizationIteration)
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

            self.__updateShape()

            self.__logCurrentOptimizationStep( optimizationIteration )

            self.__timeOptimizationStep()

            if self.__isAlgorithmConverged( optimizationIteration ):
                break        
                
    # --------------------------------------------------------------------------
    def __finalizeOptimizationLoop( self ):
        self.dataLogger.finalizeDataLogging() 

    # --------------------------------------------------------------------------
    def __callCoumminicatorToCreateNewRequests( self ):
        self.communicator.initializeCommunication()
        self.communicator.requestFunctionValueOf( self.onlyObjective )
        self.communicator.requestFunctionValueOf( self.onlyConstraint )
        self.communicator.requestGradientOf( self.onlyObjective )
        self.communicator.requestGradientOf( self.onlyConstraint )    

    # --------------------------------------------------------------------------
    def __callAnalyzerToPerformRequestedAnalyses( self, optimizationIteration ):
        self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )

    # --------------------------------------------------------------------------
    def __storeResultOfSensitivityAnalysisOnNodes( self ):
        gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( self.onlyObjective )
        gradientOfConstraintFunction = self.communicator.getReportedGradientOf ( self.onlyConstraint )
        self.__storeGradientOnNodalVariable( gradientOfObjectiveFunction, OBJECTIVE_SENSITIVITY )
        self.__storeGradientOnNodalVariable( gradientOfConstraintFunction, CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __storeGradientOnNodalVariable( self, gradients, variable_name ):
        for nodeId in gradients:
            gradient = Vector(3)
            gradient[0] = gradients[nodeId][0]
            gradient[1] = gradients[nodeId][1]
            gradient[2] = gradients[nodeId][2]           
            self.designSurface.Nodes[nodeId].SetSolutionStepValue(variable_name,0,gradient)

    # --------------------------------------------------------------------------
    def __alignSensitivitiesToLocalSurfaceNormal( self ):
            self.geometryTools.compute_unit_surface_normals()
            self.geometryTools.project_nodal_variable_on_unit_surface_normals( OBJECTIVE_SENSITIVITY )
            self.geometryTools.project_nodal_variable_on_unit_surface_normals( CONSTRAINT_SENSITIVITY )            

    # --------------------------------------------------------------------------
    def __dampSensitivities( self ):
        self.dampingUtilities.DampNodalVariable( OBJECTIVE_SENSITIVITY )
        self.dampingUtilities.DampNodalVariable( CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __computeShapeUpdate( self ):
        self.__mapSensitivitiesToDesignSpace()

        constraintValue = self.communicator.getReportedFunctionValueOf( self.onlyConstraint )
                        
        if self.__isConstraintActive( constraintValue ):
            self.optimizationTools.compute_projected_search_direction()
            self.optimizationTools.correct_projected_search_direction( constraintValue )
        else:
            self.optimizationTools.compute_search_direction_steepest_descent()

        self.optimizationTools.compute_design_update()
        self.__mapDesignUpdateToGeometrySpace() 

    # --------------------------------------------------------------------------
    def __mapSensitivitiesToDesignSpace( self ):
        self.mapper.MapToDesignSpace( OBJECTIVE_SENSITIVITY, MAPPED_OBJECTIVE_SENSITIVITY )
        self.mapper.MapToDesignSpace( CONSTRAINT_SENSITIVITY, MAPPED_CONSTRAINT_SENSITIVITY )

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
        self.mapper.MapToGeometrySpace( DESIGN_UPDATE, SHAPE_UPDATE ) 

    # --------------------------------------------------------------------------
    def __dampShapeUpdate( self ):
        self.dampingUtilities.DampNodalVariable( SHAPE_UPDATE )
        
    # --------------------------------------------------------------------------
    def __updateShape( self ):
        self.geometryTools.update_coordinates_according_to_input_variable( SHAPE_UPDATE )

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep( self, optimizationIteration ):
        self.dataLogger.logCurrentData( optimizationIteration )

    # --------------------------------------------------------------------------
    def __timeOptimizationStep( self ):
        print("\n> Time needed for current optimization step = ", self.timer.getLapTime(), "s")
        print("> Time needed for total optimization so far = ", self.timer.getTotalTime(), "s")
        self.timer.resetLapTime()

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged( self, optimizationIteration ):

        if optimizationIteration > 1 :

            # Check if maximum iterations were reached
            if optimizationIteration == self.maxIterations:
                print("\n> Maximal iterations of optimization problem reached!")
                return True

            relativeChangeOfObjectiveValue = self.dataLogger.getValue( "RELATIVE_CHANGE_OF_OBJECTIVE_VALUE" )
            
            # Check for relative tolerance
            relativeTolerance = self.optimizationSettings["optimization_algorithm"]["relative_tolerance"].GetDouble()
            if abs(relativeChangeOfObjectiveValue) < relativeTolerance:
                print("\n> Optimization problem converged within a relative objective tolerance of ",relativeTolerance,"%.")
                return True

            # Check if value of objective increases
            if relativeChangeOfObjectiveValue > 0:
                print("\n> Value of objective function increased!")
                return False

# ==============================================================================
