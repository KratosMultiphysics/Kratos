# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Long Chen https://github.com/longchentum
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.EigenSolversApplication import *

# Additional imports
from algorithm_base import OptimizationAlgorithm
import mapper_factory
import data_logger_factory
from custom_timer import Timer

# ==============================================================================
class AlgorithmSVD( OptimizationAlgorithm ) :
    # --------------------------------------------------------------------------
    def __init__( self, OptimizationSettings, ModelPartController, Analyzer, Communicator ):
        self.OptimizationSettings = OptimizationSettings
        self.ModelPartController = ModelPartController
        self.Analyzer = Analyzer
        self.Communicator = Communicator

        self.OptimizationModelPart = ModelPartController.GetOptimizationModelPart()
        self.DesignSurface = ModelPartController.GetDesignSurface()

        #self.onlyObjective = OptimizationSettings["objectives"][0]["identifier"].GetString()
        #self.onlyConstraint = OptimizationSettings["constraints"][0]["identifier"].GetString()
        self.onlyObjectiveId = OptimizationSettings["objectives"][0]["identifier"].GetString()
        self.onlyConstraintId = OptimizationSettings["constraints"][0]["identifier"].GetString()

        self.maxIterations = OptimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1
        self.projectionOnNormalsIsSpecified = OptimizationSettings["optimization_algorithm"]["project_gradients_on_surface_normals"].GetBool()

        self.dampingIsSpecified = OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()

        self.Mapper = mapper_factory.CreateMapper( ModelPartController, OptimizationSettings )
        self.DataLogger = data_logger_factory.CreateDataLogger( ModelPartController, Communicator, OptimizationSettings )

        #self.correlationFactor = OptimizationSettings["optimization_algorithm"]["correlation_factor"].GetDouble()
        self.correlationFactor = 8.0

        self.GeometryUtilities = GeometryUtilities( self.DesignSurface )
        self.OptimizationUtilities = OptimizationUtilities( self.DesignSurface, OptimizationSettings )

        self.ConvergedNodelist = []

        if self.dampingIsSpecified:
            damping_regions = self.ModelPartController.GetDampingRegions()
            self.DampingUtilities = DampingUtilities( self.DesignSurface, damping_regions, self.OptimizationSettings )

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop( self ):
        self.Analyzer.InitializeBeforeOptimizationLoop()
        self.ModelPartController.InitializeMeshController()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        timer = Timer()
        timer.StartTimer()

        for self.optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting SVD optimization iteration ",self.optimizationIteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            self.__initializeNewShape()

            self.__analyzeShape()

            constraintValue = self.Communicator.getValue( self.onlyConstraintId )

            print("###################")
            print("constraintValue:")
            print(constraintValue)
            print("###################")               

            
            if self.projectionOnNormalsIsSpecified:
                self.__projectSensitivitiesOnSurfaceNormals()

            if self.dampingIsSpecified:
                self.__dampSensitivities()

            self.__computeShapeUpdate()

            if self.dampingIsSpecified:
                self.__dampShapeUpdate()

            self.__logCurrentOptimizationStep()

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        self.DataLogger.FinalizeDataLogging()
        self.Analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewShape( self ):
        self.ModelPartController.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )
        self.ModelPartController.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf( self.onlyObjectiveId )
        self.Communicator.requestGradientOf( self.onlyObjectiveId )
        self.Communicator.requestValueOf( self.onlyConstraintId )
        self.Communicator.requestGradientOf( self.onlyConstraintId )


        self.Analyzer.AnalyzeDesignAndReportToCommunicator( self.DesignSurface, self.optimizationIteration, self.Communicator )

        self.__storeResultOfSensitivityAnalysisOnNodes()
        self.__RevertPossibleShapeModificationsDuringAnalysis()

    # # --------------------------------------------------------------------------
    # def __storeResultOfSensitivityAnalysisOnNodes( self ):
    #     gradientOfObjectiveFunction = self.Communicator.getStandardizedGradient ( self.onlyObjectiveId )
    #     gradientOfConstraintFunction = self.Communicator.getStandardizedGradient( self.onlyConstraint )
    #     for nodeId, tmp_gradient in gradientOfObjectiveFunction.items():
    #         self.OptimizationModelPart.Nodes[nodeId].SetSolutionStepValue(OBJECTIVE_SENSITIVITY,0,tmp_gradient)
    #         #self.OptimizationModelPart.N

    # --------------------------------------------------------------------------
    def __storeResultOfSensitivityAnalysisOnNodes( self ):
        self.ConvergedNodelist = self.Communicator.getReportedNodeList()
        print("__storeResultOfSensitivityAnalysisOnNodes:")
        print("ListOfConvergedNodes")
        print(self.ConvergedNodelist)
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")        
        gradientOfObjectiveFunction = self.Communicator.getStandardizedGradient ( self.onlyObjectiveId )
        gradientOfConstraintFunction = self.Communicator.getStandardizedGradient( self.onlyConstraintId )
        self.__storeGradientOnNodalVariable( gradientOfObjectiveFunction, OBJECTIVE_SENSITIVITY )
        self.__storeGradientOnNodalVariable( gradientOfConstraintFunction, CONSTRAINT_SENSITIVITY )



    # --------------------------------------------------------------------------
    def __storeGradientOnNodalVariable( self, gradients, variable_name ):
        for nodeId in gradients:
            gradient = Vector(3)
            gradient[0] = gradients[nodeId][0]
            gradient[1] = gradients[nodeId][1]
            gradient[2] = gradients[nodeId][2]
            if nodeId in self.ConvergedNodelist:
                gradient[0] = 0.0
                gradient[1] = 0.0
                gradient[2] = 0.0
                #print("TOUCH BOUND")
            #self.DesignSurface.Nodes[nodeId].SetSolutionStepValue(variable_name,0,gradient)            
            self.OptimizationModelPart.Nodes[nodeId].SetSolutionStepValue(variable_name,0,gradient)

    # --------------------------------------------------------------------------
    def __RevertPossibleShapeModificationsDuringAnalysis( self ):
        self.ModelPartController.SetMeshToReferenceMesh()
        self.ModelPartController.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def __projectSensitivitiesOnSurfaceNormals( self ):
        self.GeometryUtilities.ComputeUnitSurfaceNormals()
        self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( OBJECTIVE_SENSITIVITY )
        self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __dampSensitivities( self ):
        self.DampingUtilities.DampNodalVariable( OBJECTIVE_SENSITIVITY )
        self.DampingUtilities.DampNodalVariable( CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __computeShapeUpdate( self ):
        # self.__mapSensitivitiesToDesignSpace()
        # self.OptimizationUtilities.ComputeSearchDirectionSteepestDescent()
        # self.OptimizationUtilities.ComputeControlPointUpdate()
        # self.__mapDesignUpdateToGeometrySpace()
        
        self.__mapSensitivitiesToDesignSpace()
        constraintValue = self.Communicator.getValue( self.onlyConstraintId )


        size = self.DesignSurface.NumberOfNodes() * 3

        # print("LLLLLLLLLLLLLLLLL")
        # print(size)
        # print("LLLLLLLLLLLLLLLLL")

        M = Matrix(2, size)

        i = 0
        for node_i in self.DesignSurface.Nodes:
            dFds_i = node_i.GetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY)
            # dFds_i = node_i.GetSolutionStepValue(OBJECTIVE_SENSITIVITY)
            M[0,3*i] = dFds_i[0]
            M[0,3*i+1] = dFds_i[1]
            M[0,3*i+2] = dFds_i[2]
            dCds_i = node_i.GetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY)
            # dCds_i = node_i.GetSolutionStepValue(CONSTRAINT_SENSITIVITY)
            M[1,3*i] = dCds_i[0]
            M[1,3*i+1] = dCds_i[1]
            M[1,3*i+2] = dCds_i[2]
            #print(dCds_i)
            i = i + 1

        U = Matrix()
        V = Matrix()
        S = Matrix()

        svd_utility = SingularValueDecomposition()

        print("\n AlgorithmSVD: SVD start!!")
        svd_utility.Solve(M,U,V,S)
        print("\n AlgorithmSVD: SVD finished!!")

        a1 = 1.0
        a2 = 1.0

        if U[0,0] > 0.0:
            a1 = -1.0
        if U[0,1] > 0.0:
            a2 = -1.0

        # if self.__isConstraintActive( constraintValue ):
        #     value = self.optimizationTools.ComputeSVDSearchDirection(M,U,V, self.correlationFactor)
        #     self.DataLogger.SetValue("COS_1",value)
        # else:
        value = self.OptimizationUtilities.ComputeSVDSearchDirection(M,U,V, self.correlationFactor)
        # self.DataLogger.SetValue("COS_1",value)        

        self.OptimizationUtilities.ComputeControlPointUpdate()
        self.__mapDesignUpdateToGeometrySpace()



    # --------------------------------------------------------------------------
    def __mapSensitivitiesToDesignSpace( self ):
        self.Mapper.MapToDesignSpace( OBJECTIVE_SENSITIVITY, MAPPED_OBJECTIVE_SENSITIVITY )
        self.Mapper.MapToDesignSpace( CONSTRAINT_SENSITIVITY, MAPPED_CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __mapDesignUpdateToGeometrySpace( self ):
        self.Mapper.MapToGeometrySpace( CONTROL_POINT_UPDATE, SHAPE_UPDATE )

    # --------------------------------------------------------------------------
    def __dampShapeUpdate( self ):
        self.DampingUtilities.DampNodalVariable( SHAPE_UPDATE )

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep( self ):
        self.DataLogger.LogCurrentData( self.optimizationIteration )

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged( self ):

        if self.optimizationIteration > 1:

            # Check if maximum iterations were reached
            if self.optimizationIteration == self.maxIterations:
                print("\n> Maximal iterations of optimization problem reached!")
                return True

            relativeChangeOfObjectiveValue = self.DataLogger.GetValue( "RELATIVE_CHANGE_OF_OBJECTIVE_VALUE" )

            # Check for relative tolerance
            relativeTolerance = self.OptimizationSettings["optimization_algorithm"]["relative_tolerance"].GetDouble()
            if abs(relativeChangeOfObjectiveValue) < relativeTolerance:
                print("\n> Optimization problem converged within a relative objective tolerance of ",relativeTolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges( self ):
        self.OptimizationUtilities.AddFirstVariableToSecondVariable( CONTROL_POINT_UPDATE, CONTROL_POINT_CHANGE )
        self.OptimizationUtilities.AddFirstVariableToSecondVariable( SHAPE_UPDATE, SHAPE_CHANGE )            

# ==============================================================================




