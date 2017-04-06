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
        self.correctionScaling = optimizationSettings["optimization_algorithm"]["correction_scaling"].GetDouble()
        self.stepSize = optimizationSettings["line_search"]["step_size"].GetDouble()


        self.geometryTools = GeometryUtilities( designSurface )
        self.optimizationTools = OptimizationUtilities( designSurface, optimizationSettings )
        self.timer = (__import__("timer_factory")).CreateTimer()

        additionalVariablesToBeLogged = { "stepSize": self.stepSize, "correctionScaling": self.correctionScaling }
        self.dataLogger = (__import__("optimization_data_logger_factory")).CreateDataLogger( designSurface, 
                                                                                             communicator, 
                                                                                             optimizationSettings, 
                                                                                             self.timer, 
                                                                                             additionalVariablesToBeLogged  )        
    # --------------------------------------------------------------------------
    def executeAlgorithm( self ):
        self.initializeOptimizationLoop()
        self.startOptimizationLoop()
        self.finalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def initializeOptimizationLoop( self ):   
        self.timer.startTimer()
        self.dataLogger.initializeDataLogging() 
        self.optimizationTools.set_correction_scaling( self.correctionScaling )

    # --------------------------------------------------------------------------
    def startOptimizationLoop( self ):

        for optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",self.timer.getTimeStamp(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            self.communicator.initializeCommunication()
            self.communicator.requestFunctionValueOf( self.onlyObjective )
            self.communicator.requestFunctionValueOf( self.onlyConstraint )
            self.communicator.requestGradientOf( self.onlyObjective )
            self.communicator.requestGradientOf( self.onlyConstraint )

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )

            constraintValue = self.communicator.getReportedFunctionValueOf ( self.onlyConstraint )
            constraintIsActive = self.checkIfConstraintIsActive( constraintValue )
            self.storeGradientObtainedByCommunicatorOnNodes()

            self.geometryTools.compute_unit_surface_normals()
            self.geometryTools.project_grad_on_unit_surface_normal( constraintIsActive )

            self.mapper.compute_mapping_matrix()
            self.mapper.map_sensitivities_to_design_space( constraintIsActive )

            if constraintIsActive:
                self.optimizationTools.compute_projected_search_direction( constraintValue )
                self.optimizationTools.correct_projected_search_direction( constraintValue )
                self.correctionScaling = self.optimizationTools.get_correction_scaling()
            else:
                self.optimizationTools.compute_search_direction_steepest_descent()
                self.correctionScaling = "-"

            self.optimizationTools.compute_design_update()
            self.mapper.map_design_update_to_geometry_space()  

            self.dataLogger.logCurrentData( optimizationIteration )

            print("\n> Time needed for current optimization step = ", self.timer.getLapTime(), "s")
            print("> Time needed for total optimization so far = ", self.timer.getTotalTime(), "s") 
            self.timer.resetLapTime()

            if self.algorithmConverged( optimizationIteration ):
                break

    # --------------------------------------------------------------------------
    def finalizeOptimizationLoop( self ):
        self.dataLogger.finalizeDataLogging() 

    # --------------------------------------------------------------------------
    def storeGradientObtainedByCommunicatorOnNodes( self ):

        gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( self.onlyObjective )
        gradientOfConstraintFunction = {}
        if self.onlyConstraint != None:
            gradientOfConstraintFunction = self.communicator.getReportedGradientOf ( self.onlyConstraint )

        # Read objective gradients
        for node_Id in gradientOfObjectiveFunction:

            # If deactivated, nodal sensitivities will not be assigned and hence remain zero
            if self.designSurface.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                continue

            # If not deactivated, nodal sensitivities will be assigned
            sens_i = Vector(3)
            sens_i[0] = gradientOfObjectiveFunction[node_Id][0]
            sens_i[1] = gradientOfObjectiveFunction[node_Id][1]
            sens_i[2] = gradientOfObjectiveFunction[node_Id][2]           
            self.designSurface.Nodes[node_Id].SetSolutionStepValue(OBJECTIVE_SENSITIVITY,0,sens_i)

        # When gradientOfConstraintFunction is defined also store constraint sensitivities (bool returns false if dictionary is empty)
        if self.onlyConstraint != None:
            eucledian_norm_cons_sens = 0.0
            for node_Id in gradientOfConstraintFunction:

                # If deactivated, nodal sensitivities will not be assigned and hence remain zero
                if self.designSurface.Nodes[node_Id].GetSolutionStepValue(SENSITIVITIES_DEACTIVATED):
                    continue

                # If not deactivated, nodal sensitivities will be assigned
                sens_i = Vector(3)
                sens_i[0] = gradientOfConstraintFunction[node_Id][0]
                sens_i[1] = gradientOfConstraintFunction[node_Id][1]
                sens_i[2] = gradientOfConstraintFunction[node_Id][2]           
                self.designSurface.Nodes[node_Id].SetSolutionStepValue(CONSTRAINT_SENSITIVITY,0,sens_i)     

    # --------------------------------------------------------------------------
    def checkIfConstraintIsActive( self, constraintValue ):
        if self.typOfOnlyConstraint == "equality":
            return True
        elif constraintValue > 0:
            return True
        else:
            return False              

    # --------------------------------------------------------------------------
    def algorithmConverged( self, optimizationIteration ):

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
