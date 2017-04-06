# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
def CreateOptimizer( inputModelPart, optimizationSettings ):
    if  optimizationSettings["design_variables"]["variable_type"].GetString() == "vertex_morphing":
        optimizer = VertexMorphingMethod( inputModelPart, optimizationSettings )
        return optimizer
    else:
        sys.exit( "Specified design_control not implemented" )

# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__( self, inputModelPart, optimizationSettings ):
        
        self.inputModelPart = inputModelPart
        self.optimizationSettings = optimizationSettings
        self.addVariablesNeededForOptimization( inputModelPart )

    # --------------------------------------------------------------------------
    def addVariablesNeededForOptimization( self, inputModelPart ):
        inputModelPart.AddNodalSolutionStepVariable(NORMAL)
        inputModelPart.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY) 
        inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY) 
        inputModelPart.AddNodalSolutionStepVariable(DESIGN_UPDATE)
        inputModelPart.AddNodalSolutionStepVariable(DESIGN_CHANGE_ABSOLUTE)  
        inputModelPart.AddNodalSolutionStepVariable(SEARCH_DIRECTION) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATE) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        inputModelPart.AddNodalSolutionStepVariable(IS_ON_BOUNDARY)
        inputModelPart.AddNodalSolutionStepVariable(BOUNDARY_PLANE) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATES_DEACTIVATED) 
        inputModelPart.AddNodalSolutionStepVariable(SENSITIVITIES_DEACTIVATED) 

    # --------------------------------------------------------------------------
    def importModelPart( self ):
        model_part_io = ModelPartIO( self.optimizationSettings["design_variables"]["input_model_part_name"].GetString() )
        model_part_io.ReadModelPart( self.inputModelPart )
        buffer_size = 1
        self.inputModelPart.SetBufferSize( buffer_size )
        self.inputModelPart.ProcessInfo.SetValue( DOMAIN_SIZE, self.optimizationSettings["design_variables"]["domain_size"].GetInt() )

    # --------------------------------------------------------------------------
    def importAnalyzer( self, newAnalyzer ): 
        self.analyzer = newAnalyzer

    # --------------------------------------------------------------------------
    def optimize( self ):
        
        # timer
        timerFactory = __import__("timer_factory")
        self.timer = timerFactory.CreateTimer()

        print("\n> ==============================================================================================================")
        print("> ",self.timer.getTimeStamp(),": Starting optimization using the following algorithm: ",self.optimizationSettings["optimization_algorithm"]["name"].GetString())
        print("> ==============================================================================================================\n")
    
        # model data
        self.designSurface = self.getDesignSurfaceFromInputModelPart()
        self.listOfDampingRegions = self.getListOfDampingRegionsFromInputModelPart()

        # mapper
        self.vertexMorphingMapper = VertexMorphingMapper( self.designSurface, self.listOfDampingRegions, self.optimizationSettings["design_variables"] ) 

        # communicator
        communicatorFactory = __import__("communicator_factory")        
        self.communicator = communicatorFactory.CreateCommunicator( self.optimizationSettings )    

        # algorithm_driver
        # .....algorithmDriver(timer, communicator, design_surface, vertexMorphingMapper, optimizationDataLogger )

        # data logger
        optimizationDataLoggerFactory = __import__("optimization_data_logger_factory")
        self.optimizationDataLogger = optimizationDataLoggerFactory.CreateDataLogger( self.designSurface, self.communicator, self.timer, self.optimizationSettings )

        self.optimizationDataLogger.initializeDataLogging() 
        self.runSpecifiedOptimizationAlgorithm()
        self.optimizationDataLogger.finalizeDataLogging()        

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")
    
    # --------------------------------------------------------------------------
    def getDesignSurfaceFromInputModelPart( self ):
        nameOfDesingSurface = self.optimizationSettings["design_variables"]["design_submodel_part_name"].GetString()
        if self.inputModelPart.HasSubModelPart( nameOfDesingSurface ):
            optimizationModel = self.inputModelPart.GetSubModelPart( nameOfDesingSurface )
            print("> The following design surface was defined:\n\n",optimizationModel)
            return optimizationModel
        else:
            raise RuntimeError("Sub-model part specified for optimization does not exist!")

    # --------------------------------------------------------------------------
    def getListOfDampingRegionsFromInputModelPart( self ):
        listOfDampingRegions = {}
        print("> The following damping regions are defined: \n")
        for regionNumber in range(self.optimizationSettings["design_variables"]["damping"]["damping_regions"].size()):
            regionName = self.optimizationSettings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
            if self.inputModelPart.HasSubModelPart(regionName):
                print(regionName)
                listOfDampingRegions[regionName] = self.inputModelPart.GetSubModelPart(regionName)
            else:
                raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)         
        return listOfDampingRegions

    # --------------------------------------------------------------------------
    def runSpecifiedOptimizationAlgorithm( self ):
        if self.optimizationSettings["optimization_algorithm"]["name"].GetString() == "steepest_descent":
           self.runSteepestDescentAlgorithm()
        elif self.optimizationSettings["optimization_algorithm"]["name"].GetString() == "penalized_projection":
           self.runPenalizedProjectionAlgorithm()           
        else:
            sys.exit("Specified optimization_algorithm not implemented!")

    # --------------------------------------------------------------------------
    def runSteepestDescentAlgorithm( self ):

        # Flags to trigger proper function calls
        constraints_given = False

        # Get Id of objective (right now only one objective is considered!!)
        onlyObjective = self.optimizationSettings["objectives"][0]["identifier"].GetString()

        # Tools run optimization algorithm
        geometryTools = GeometryUtilities( self.designSurface )
        optimizationTools = OptimizationUtilities( self.designSurface, self.optimizationSettings )

        # Start optimization loop
        self.timer.startTimer()
        maxOptimizationIterations = self.optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1
        for optimizationIteration in range(1, maxOptimizationIterations):

            # Some output
            print("\n>===================================================================")
            print("> ",self.timer.getTimeStamp(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            self.communicator.initializeCommunication()
            self.communicator.requestFunctionValueOf( onlyObjective )
            self.communicator.requestGradientOf( onlyObjective )

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )
         
            self.storeGradientObtainedByCommunicatorOnNodes( onlyObjective )

            geometryTools.compute_unit_surface_normals()
            geometryTools.project_grad_on_unit_surface_normal( constraints_given )

            self.vertexMorphingMapper.compute_mapping_matrix()
            self.vertexMorphingMapper.map_sensitivities_to_design_space( constraints_given )

            optimizationTools.compute_search_direction_steepest_descent()
            optimizationTools.compute_design_update()

            self.vertexMorphingMapper.map_design_update_to_geometry_space()         
            
            self.optimizationDataLogger.logCurrentData( optimizationIteration )

            # Take time needed for current optimization step
            print("\n> Time needed for current optimization step = ", self.timer.getLapTime(), "s")
            print("> Time needed for total optimization so far = ", self.timer.getTotalTime(), "s") 

            # Check convergence
            if optimizationIteration > 1 :

                # Check if maximum iterations were reached
                if optimizationIteration == maxOptimizationIterations:
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                relativeChangeOfObjectiveValue = self.optimizationDataLogger.getRelativeChangeOfObjectiveValue( optimizationIteration )
                
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
    def runPenalizedProjectionAlgorithm( self ):

        # Flags to trigger proper function calls
        constraints_given = True

        # Get information about response functions
        onlyObjective = self.optimizationSettings["objectives"][0]["identifier"].GetString()
        onlyConstraint = self.optimizationSettings["constraints"][0]["identifier"].GetString()
        typOfOnlyConstraint = self.optimizationSettings["constraints"][0]["type"].GetString()      

        # Tools run optimization algorithm
        geometryTools = GeometryUtilities( self.designSurface )
        optimizationTools = OptimizationUtilities( self.designSurface, self.optimizationSettings )      

        # Start optimization loop
        self.timer.startTimer()
        maxOptimizationIterations = self.optimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1
        for optimizationIteration in range(1,maxOptimizationIterations):

            # Some output
            print("\n>===================================================================")
            print("> ",self.timer.getTimeStamp(),": Starting optimization iteration ",optimizationIteration)
            print(">===================================================================\n")

            self.communicator.initializeCommunication()
            self.communicator.requestFunctionValueOf( onlyObjective )
            self.communicator.requestFunctionValueOf( onlyConstraint )
            self.communicator.requestGradientOf( onlyObjective )
            self.communicator.requestGradientOf( onlyConstraint )

            self.analyzer.analyzeDesignAndReportToCommunicator( self.designSurface, optimizationIteration, self.communicator )

            constraintValue = self.communicator.getReportedFunctionValueOf ( onlyConstraint )

            self.storeGradientObtainedByCommunicatorOnNodes( onlyObjective, onlyConstraint )
            
            # Check if constraint is active
            if typOfOnlyConstraint == "equality":
                constraints_given = True
            elif constraintValue > 0:
                constraints_given = True
            else:
                constraints_given = False       

            geometryTools.compute_unit_surface_normals()
            geometryTools.project_grad_on_unit_surface_normal( constraints_given )

            self.vertexMorphingMapper.compute_mapping_matrix()
            self.vertexMorphingMapper.map_sensitivities_to_design_space( constraints_given )

            correctionScaling = [False] 
            if constraints_given:
                optimizationTools.compute_projected_search_direction( constraintValue )
                optimizationTools.correct_projected_search_direction( constraintValue, correctionScaling )
                optimizationTools.compute_design_update()
            else:
                optimizationTools.compute_search_direction_steepest_descent()
                optimizationTools.compute_design_update()

            self.vertexMorphingMapper.map_design_update_to_geometry_space()  

            self.optimizationDataLogger.logCurrentData( optimizationIteration )

            # Take time needed for current optimization step
            print("\n> Time needed for current optimization step = ", self.timer.getLapTime(), "s")
            print("> Time needed for total optimization so far = ", self.timer.getTotalTime(), "s") 

            # Check convergence
            if optimizationIteration > 1 :

                # Check if maximum iterations were reached
                if optimizationIteration == maxOptimizationIterations:
                    print("\n> Maximal iterations of optimization problem reached!")
                    break

                relativeChangeOfObjectiveValue = self.optimizationDataLogger.getRelativeChangeOfObjectiveValue( optimizationIteration )
                
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
    def storeGradientObtainedByCommunicatorOnNodes( self, onlyObjective, onlyConstraint=None ):

        gradientOfObjectiveFunction = self.communicator.getReportedGradientOf ( onlyObjective )
        gradientOfConstraintFunction = {}
        if onlyConstraint != None:
            gradientOfConstraintFunction = self.communicator.getReportedGradientOf ( onlyConstraint )

        # Read objective gradients
        eucledian_norm_obj_sens = 0.0
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
        if onlyConstraint != None:
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

# ==============================================================================