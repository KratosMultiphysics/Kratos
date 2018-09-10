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

# Additional imports
import shutil
import os

from design_logger_gid import DesignLoggerGID
from design_logger_unv import DesignLoggerUNV
from design_logger_vtk import DesignLoggerVTK

from value_logger_steepest_descent import ValueLoggerSteepestDescent
from value_logger_penalized_projection import ValueLoggerPenalizedProjection
from value_logger_trust_region import ValueLoggerTrustRegion

# ==============================================================================
def CreateDataLogger( ModelPartController, Communicator, OptimizationSettings ):
    return DataLogger( ModelPartController, Communicator, OptimizationSettings )

# ==============================================================================
class DataLogger():
    # --------------------------------------------------------------------------
    def __init__( self, ModelPartController, Communicator, OptimizationSettings ):
        self.ModelPartController = ModelPartController
        self.Communicator = Communicator
        self.OptimizationSettings = OptimizationSettings

        self.ValueLogger = self.__CreateValueLogger()
        self.DesignLogger = self.__CreateDesignLogger()

        self.__CreateFolderToStoreOptimizationResults()
        self.__OutputInformationAboutResponseFunctions()

    # -----------------------------------------------------------------------------
    def __CreateValueLogger( self ):
        AlgorithmName = self.OptimizationSettings["optimization_algorithm"]["name"].GetString()
        if AlgorithmName == "steepest_descent":
            return ValueLoggerSteepestDescent( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "penalized_projection":
            return ValueLoggerPenalizedProjection( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "trust_region":
            return ValueLoggerTrustRegion( self.Communicator, self.OptimizationSettings )
        else:
            raise NameError("The following optimization algorithm not supported by the response logger (name may be misspelled): " + AlgorithmName)

    # -----------------------------------------------------------------------------
    def __CreateDesignLogger( self ):
        outputFormatName = self.OptimizationSettings["output"]["output_format"]["name"].GetString()
        if outputFormatName == "gid":
            return DesignLoggerGID( self.ModelPartController, self.OptimizationSettings )
        if outputFormatName == "unv":
            return DesignLoggerUNV( self.ModelPartController, self.OptimizationSettings )
        if outputFormatName == "vtk":
            return DesignLoggerVTK( self.ModelPartController, self.OptimizationSettings )
        else:
            raise NameError("The following output format is not supported by the design logger (name may be misspelled): " + outputFormatName)

    # --------------------------------------------------------------------------
    def __CreateFolderToStoreOptimizationResults ( self ):
        resultsDirectory = self.OptimizationSettings["output"]["output_directory"].GetString()
        if os.path.exists(resultsDirectory):
            shutil.rmtree(resultsDirectory)
        os.makedirs(resultsDirectory)

    # --------------------------------------------------------------------------
    def __OutputInformationAboutResponseFunctions( self ):
        numberOfObjectives = self.OptimizationSettings["objectives"].size()
        numberOfConstraints = self.OptimizationSettings["constraints"].size()

        print("\n> The following objectives are defined:\n")
        for objectiveNumber in range(numberOfObjectives):
            print(self.OptimizationSettings["objectives"][objectiveNumber])

        if numberOfConstraints != 0:
            print("> The following constraints are defined:\n")
            for constraintNumber in range(numberOfConstraints):
                print(self.OptimizationSettings["constraints"][constraintNumber],"\n")
        else:
            print("> No constraints defined.\n")

    # --------------------------------------------------------------------------
    def InitializeDataLogging( self ):
        self.DesignLogger.InitializeLogging()
        self.ValueLogger.InitializeLogging()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, current_iteration ):
        self.DesignLogger.LogCurrentDesign( current_iteration )

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, current_iteration, additional_values ):
        self.ValueLogger.LogCurrentValues( current_iteration, additional_values )

    # --------------------------------------------------------------------------
    def FinalizeDataLogging( self ):
        self.DesignLogger.FinalizeLogging()
        self.ValueLogger.FinalizeLogging()

    # --------------------------------------------------------------------------
    def GetValue( self, key, iteration ):
        return self.ValueLogger.GetValue(key, iteration)

    # --------------------------------------------------------------------------
    def GetValueHistory( self, key ):
        return self.ValueLogger.GetValueHistory(key)

# ==============================================================================
