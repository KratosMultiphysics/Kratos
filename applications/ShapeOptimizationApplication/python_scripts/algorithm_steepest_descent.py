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
    def __init__( self, designSurface, analyzer, mapper, communicator, dataLogger, timer, optimizationSettings ):

        self.designSurface = designSurface
        self.analyzer = analyzer
        self.mapper = mapper
        self.communicator = communicator
        self.dataLogger = dataLogger
        self.timer = timer 
        self.optimizationSettings = optimizationSettings

        self.maxIterations = optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1        
        self.onlyObjective = optimizationSettings["objectives"][0]["identifier"].GetString()
        self.constraints_given = False

        self.geometryTools = GeometryUtilities( designSurface )
        self.optimizationTools = OptimizationUtilities( designSurface, optimizationSettings )

    # --------------------------------------------------------------------------
    def executeAlgorithm( self ):
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

            self.communicator.initializeCommunication()
            self.communicator.requestFunctionValueOf( self.onlyObjective )
            self.communicator.requestGradientOf( self.onlyObjective )

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )
         
            self.storeGradientObtainedByCommunicatorOnNodes()

            self.geometryTools.compute_unit_surface_normals()
            self.geometryTools.project_grad_on_unit_surface_normal( self.constraints_given )

            self.mapper.compute_mapping_matrix()
            self.mapper.map_sensitivities_to_design_space( self.constraints_given )

            self.optimizationTools.compute_search_direction_steepest_descent()
            self.optimizationTools.compute_design_update()

            self.mapper.map_design_update_to_geometry_space()         
            
            self.dataLogger.logCurrentData( optimizationIteration )

            # Take time needed for current optimization step
            print("\n> Time needed for current optimization step = ", self.timer.getLapTime(), "s")
            print("> Time needed for total optimization so far = ", self.timer.getTotalTime(), "s") 

            # Check convergence
            if optimizationIteration > 1 :

                # Check if maximum iterations were reached
                if optimizationIteration == self.maxIterations:
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                relativeChangeOfObjectiveValue = self.dataLogger.getRelativeChangeOfObjectiveValue( optimizationIteration )
                
                # Check for relative tolerance
                relativeTolerance = self.optimizationSettings["optimization_algorithm"]["relative_tolerance"].GetDouble()
                if abs(relativeChangeOfObjectiveValue) < relativeTolerance:
                    print("\n> Optimization problem converged within a relative objective tolerance of ",relativeTolerance,"%.")
                    break

                # Check if value of objective increases
                if relativeChangeOfObjectiveValue > 0:
                    print("\n> Value of objective function increased!")

            self.timer.resetLapTime()

    # --------------------------------------------------------------------------
    def finalizeOptimizationLoop( self ):
        self.dataLogger.finalizeDataLogging() 

    # --------------------------------------------------------------------------
    def storeGradientObtainedByCommunicatorOnNodes( self ):

        gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( self.onlyObjective )

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
                        
# ==============================================================================
