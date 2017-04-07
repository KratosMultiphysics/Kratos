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

# ==============================================================================
class AlgorithmSteepestDescent( OptimizationAlgorithm ) :

    # --------------------------------------------------------------------------
    def __init__( self, designSurface, analyzer, mapper, communicator, optimizationSettings ):

        self.designSurface = designSurface
        self.analyzer = analyzer
        self.mapper = mapper
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings

        self.maxIterations = optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1        
        self.onlyObjective = optimizationSettings["objectives"][0]["identifier"].GetString()
        self.stepSize = optimizationSettings["line_search"]["step_size"].GetDouble()
        self.isConstraintGiven = False

        self.geometryTools = GeometryUtilities( designSurface )
        self.optimizationTools = OptimizationUtilities( designSurface, optimizationSettings )

        self.timer = (__import__("timer_factory")).CreateTimer()
        specificVariablesToBeLogged = { "stepSize": self.stepSize }
        self.dataLogger = (__import__("optimization_data_logger_factory")).CreateDataLogger( designSurface, 
                                                                                             communicator, 
                                                                                             optimizationSettings, 
                                                                                             self.timer, 
                                                                                             specificVariablesToBeLogged  )             

    # --------------------------------------------------------------------------
    def execute( self ):
        self.initializeOptimizationLoop()
        self.startOptimizationLoop()
        self.finalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def initializeOptimizationLoop( self ):   
        self.timer.startTimer()
        self.dataLogger.initializeDataLogging() 

    # --------------------------------------------------------------------------
    def startOptimizationLoop( self ):

        for optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",self.timer.getTimeStamp(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            self.callCoumminicatorToCreateNewRequests()

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )
         
            self.storeGradientsObtainedByCommunicatorOnNodes()

            self.geometryTools.compute_unit_surface_normals()
            self.geometryTools.project_grad_on_unit_surface_normal( self.isConstraintGiven )

            self.mapper.compute_mapping_matrix()
            self.mapper.map_sensitivities_to_design_space( self.isConstraintGiven )

            self.optimizationTools.compute_search_direction_steepest_descent()

            self.optimizationTools.compute_design_update()

            self.mapper.map_design_update_to_geometry_space()         
            
            self.dataLogger.logCurrentData( optimizationIteration )

            # Take time needed for current optimization step
            print("\n> Time needed for current optimization step = ", self.timer.getLapTime(), "s")
            print("> Time needed for total optimization so far = ", self.timer.getTotalTime(), "s")
            self.timer.resetLapTime()

            if self.isAlgorithmConverged( optimizationIteration ):
                break            

    # --------------------------------------------------------------------------
    def finalizeOptimizationLoop( self ):
        self.dataLogger.finalizeDataLogging() 

    # --------------------------------------------------------------------------
    def callCoumminicatorToCreateNewRequests( self ):
        self.communicator.initializeCommunication()
        self.communicator.requestFunctionValueOf( self.onlyObjective )
        self.communicator.requestGradientOf( self.onlyObjective )    

    # --------------------------------------------------------------------------
    def storeGradientsObtainedByCommunicatorOnNodes( self ):
        gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( self.onlyObjective )        
        for nodeId in gradientOfObjectiveFunction:
            if self.designSurface.Nodes[nodeId].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                continue
            sensitivity = Vector(3)
            sensitivity[0] = gradientOfObjectiveFunction[nodeId][0]
            sensitivity[1] = gradientOfObjectiveFunction[nodeId][1]
            sensitivity[2] = gradientOfObjectiveFunction[nodeId][2]           
            self.designSurface.Nodes[nodeId].SetSolutionStepValue(OBJECTIVE_SENSITIVITY,0,sensitivity)

    # --------------------------------------------------------------------------
    def isAlgorithmConverged( self, optimizationIteration ):

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
