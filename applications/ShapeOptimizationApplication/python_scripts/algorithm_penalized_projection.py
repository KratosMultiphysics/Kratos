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
class AlgorithmPenalizedProjection( OptimizationAlgorithm ) :

    # --------------------------------------------------------------------------
    def __init__( self, designSurface, analyzer, mapper, communicator, optimizationSettings ):

        self.designSurface = designSurface
        self.analyzer = analyzer
        self.mapper = mapper
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings

        self.onlyObjective = optimizationSettings["objectives"][0]["identifier"].GetString()
        self.onlyConstraint = optimizationSettings["constraints"][0]["identifier"].GetString()
        self.typOfOnlyConstraint = optimizationSettings["constraints"][0]["type"].GetString()          
        self.maxIterations = optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1        
        self.initialCorrectionScaling = optimizationSettings["optimization_algorithm"]["correction_scaling"].GetDouble()
        self.initialStepSize = optimizationSettings["line_search"]["step_size"].GetDouble()

        self.geometryTools = GeometryUtilities( designSurface )
        self.optimizationTools = OptimizationUtilities( designSurface, optimizationSettings )

        self.timer = (__import__("timer_factory")).CreateTimer()
        self.specificVariablesToBeLogged = { "correctionScaling": self.initialCorrectionScaling, "stepSize": self.initialStepSize }
        self.dataLogger = (__import__("optimization_data_logger_factory")).CreateDataLogger( designSurface, 
                                                                                             communicator, 
                                                                                             optimizationSettings, 
                                                                                             self.timer, 
                                                                                             self.specificVariablesToBeLogged  )        
    # --------------------------------------------------------------------------
    def execute( self ):
        self.initializeOptimizationLoop()
        self.startOptimizationLoop()
        self.finalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def initializeOptimizationLoop( self ):   
        self.timer.startTimer()
        self.dataLogger.initializeDataLogging() 
        self.optimizationTools.set_correction_scaling( self.initialCorrectionScaling )

    # --------------------------------------------------------------------------
    def startOptimizationLoop( self ):

        for optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",self.timer.getTimeStamp(),": Starting optimization iteration ", optimizationIteration)
            print(">===================================================================\n")

            self.callCoumminicatorToCreateNewRequests()

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )

            constraintIsActive = self.isConstraintActive()
            self.storeGradientsObtainedByCommunicatorOnNodes()

            self.geometryTools.compute_unit_surface_normals()
            self.geometryTools.project_grad_on_unit_surface_normal( constraintIsActive )

            self.mapper.compute_mapping_matrix()
            self.mapper.map_sensitivities_to_design_space( constraintIsActive )

            self.computeSearchDirection( constraintIsActive )

            self.optimizationTools.compute_design_update()

            self.mapper.map_design_update_to_geometry_space()  

            self.updateSpecificVariablesForLog()            
            self.dataLogger.logCurrentData( optimizationIteration )


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
        self.communicator.requestFunctionValueOf( self.onlyConstraint )
        self.communicator.requestGradientOf( self.onlyObjective )
        self.communicator.requestGradientOf( self.onlyConstraint )    

    # --------------------------------------------------------------------------
    def storeGradientsObtainedByCommunicatorOnNodes( self ):
        self.storeObjectiveGradients()
        if self.onlyConstraint != None:
            self.storeConstraintGradients()

    # --------------------------------------------------------------------------
    def storeObjectiveGradients( self ):
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
    def storeConstraintGradients( self ):
        gradientOfConstraintFunction = self.communicator.getReportedGradientOf ( self.onlyConstraint )        
        for nodeId in gradientOfConstraintFunction:
            if self.designSurface.Nodes[nodeId].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                continue
            sensitivity = Vector(3)
            sensitivity[0] = gradientOfConstraintFunction[nodeId][0]
            sensitivity[1] = gradientOfConstraintFunction[nodeId][1]
            sensitivity[2] = gradientOfConstraintFunction[nodeId][2]           
            self.designSurface.Nodes[nodeId].SetSolutionStepValue(CONSTRAINT_SENSITIVITY,0,sensitivity)     

    # --------------------------------------------------------------------------
    def isConstraintActive( self ):
        constraintValue = self.communicator.getReportedFunctionValueOf ( self.onlyConstraint )        
        if self.typOfOnlyConstraint == "equality":
            return True
        elif constraintValue > 0:
            return True
        else:
            return False              

    # --------------------------------------------------------------------------
    def computeSearchDirection( self, constraintIsActive ):
        constraintValue = self.communicator.getReportedFunctionValueOf ( self.onlyConstraint )        
        if constraintIsActive:
            self.optimizationTools.compute_projected_search_direction( constraintValue )
            self.optimizationTools.correct_projected_search_direction( constraintValue )
        else:
            self.optimizationTools.compute_search_direction_steepest_descent()

    # --------------------------------------------------------------------------
    def updateSpecificVariablesForLog( self ):
        self.specificVariablesToBeLogged["correctionScaling"] = self.optimizationTools.get_correction_scaling()

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
