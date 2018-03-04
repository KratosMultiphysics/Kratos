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

# ==============================================================================
class AlgorithmPenalizedProjection( OptimizationAlgorithm ) :

    # --------------------------------------------------------------------------
    def __init__( self,
                  ModelPartController,
                  Analyzer,
                  Communicator,
                  Mapper,
                  DataLogger,
                  OptimizationSettings ):

        self.ModelPartController = ModelPartController
        self.Analyzer = Analyzer
        self.Communicator = Communicator
        self.Mapper = Mapper
        self.DataLogger = DataLogger
        self.OptimizationSettings = OptimizationSettings

        self.OptimizationModelPart = ModelPartController.GetOptimizationModelPart()
        self.DesignSurface = ModelPartController.GetDesignSurface()

        self.onlyObjectiveId = OptimizationSettings["objectives"][0]["identifier"].GetString()
        self.onlyConstraintId = OptimizationSettings["constraints"][0]["identifier"].GetString()
        self.typeOfOnlyConstraint = OptimizationSettings["constraints"][0]["type"].GetString()
        self.projectionOnNormalsIsSpecified = OptimizationSettings["optimization_algorithm"]["project_gradients_on_surface_normals"].GetBool()
        self.dampingIsSpecified = OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()
        self.maxIterations = OptimizationSettings["optimization_algorithm"]["max_iterations"].GetInt() + 1

        self.GeometryUtilities = GeometryUtilities( self.DesignSurface )
        self.OptimizationUtilities = OptimizationUtilities( self.DesignSurface, OptimizationSettings )
        if self.dampingIsSpecified:
            damping_regions = self.ModelPartController.GetDampingRegions()
            self.DampingUtilities = DampingUtilities( self.DesignSurface, damping_regions, self.OptimizationSettings )

    # --------------------------------------------------------------------------
    def execute( self ):
        self.__initializeOptimizationLoop()
        self.__runOptimizationLoop()
        self.__finalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeOptimizationLoop( self ):
        self.Analyzer.initializeBeforeOptimizationLoop()
        self.ModelPartController.InitializeMeshController()
        self.DataLogger.StartTimer()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def __runOptimizationLoop( self ):

        for self.optimizationIteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",self.DataLogger.GetTimeStamp(),": Starting optimization iteration ", self.optimizationIteration)
            print(">===================================================================\n")

            self.__initializeNewShape()

            self.__analyzeShape()

            if self.projectionOnNormalsIsSpecified:
                self.__projectSensitivitiesOnSurfaceNormals()

            if self.dampingIsSpecified:
                self.__dampSensitivities()

            self.__computeShapeUpdate()

            if self.dampingIsSpecified:
                self.__dampShapeUpdate()

            self.__logCurrentOptimizationStep()

            self.__timeOptimizationStep()

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def __finalizeOptimizationLoop( self ):
        self.DataLogger.FinalizeDataLogging()
        self.Analyzer.finalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewShape( self ):
        self.ModelPartController.CloneTimeStep( self.optimizationIteration )
        self.ModelPartController.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )
        self.ModelPartController.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf( self.onlyObjectiveId )
        self.Communicator.requestValueOf( self.onlyConstraintId )
        self.Communicator.requestGradientOf( self.onlyObjectiveId )
        self.Communicator.requestGradientOf( self.onlyConstraintId )

        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.optimizationIteration, self.Communicator )

        self.__storeResultOfSensitivityAnalysisOnNodes()
        self.__RevertPossibleShapeModificationsDuringAnalysis()

    # --------------------------------------------------------------------------
    def __storeResultOfSensitivityAnalysisOnNodes( self ):
        gradientOfObjectiveFunction = self.Communicator.getStandardizedGradient( self.onlyObjectiveId )
        gradientOfConstraintFunction = self.Communicator.getStandardizedGradient( self.onlyConstraintId )
        self.__storeGradientOnNodalVariable( gradientOfObjectiveFunction, OBJECTIVE_SENSITIVITY )
        self.__storeGradientOnNodalVariable( gradientOfConstraintFunction, CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __storeGradientOnNodalVariable( self, gradient, variable_name ):
        for nodeId, tmp_gradient in gradient.items():
            self.OptimizationModelPart.Nodes[nodeId].SetSolutionStepValue(variable_name,0,tmp_gradient)

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
        self.__mapSensitivitiesToDesignSpace()
        constraint_value = self.Communicator.getStandardizedValue( self.onlyConstraintId )
        if self.__isConstraintActive( constraint_value ):
            self.OptimizationUtilities.ComputeProjectedSearchDirection()
            self.OptimizationUtilities.CorrectProjectedSearchDirection( constraint_value )
        else:
            self.OptimizationUtilities.ComputeSearchDirectionSteepestDescent()
        self.OptimizationUtilities.ComputeControlPointUpdate()
        self.__mapDesignUpdateToGeometrySpace()

    # --------------------------------------------------------------------------
    def __mapSensitivitiesToDesignSpace( self ):
        self.Mapper.MapToDesignSpace( OBJECTIVE_SENSITIVITY, MAPPED_OBJECTIVE_SENSITIVITY )
        self.Mapper.MapToDesignSpace( CONSTRAINT_SENSITIVITY, MAPPED_CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __isConstraintActive( self, constraintValue ):
        if self.typeOfOnlyConstraint == "=":
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
        self.DampingUtilities.DampNodalVariable( SHAPE_UPDATE )

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep( self ):
        self.DataLogger.LogCurrentData( self.optimizationIteration )

    # --------------------------------------------------------------------------
    def __timeOptimizationStep( self ):
        print("\n> Time needed for current optimization step = ", self.DataLogger.GetLapTime(), "s")
        print("> Time needed for total optimization so far = ", self.DataLogger.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged( self ):

        if self.optimizationIteration > 1 :

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
