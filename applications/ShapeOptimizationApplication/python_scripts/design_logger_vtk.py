# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#                   Sferza Massimo, https://github.com/IIIaxS
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Import logger base classes
from design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerVTK( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, ModelPartController, OptimizationSettings ):
        self.OptimizationModelPart = ModelPartController.GetOptimizationModelPart()
        self.DesignSurface = ModelPartController.GetDesignSurface()
        self.OutputSettings = OptimizationSettings["output"]

        self.__DetermineOutputMode()
        self.__CreateVTKIO()

    # --------------------------------------------------------------------------
    def __DetermineOutputMode( self ):
        OutputMode = self.OutputSettings["design_output_mode"].GetString()

        self.WriteDesignSurface = False
        self.WriteOptimizationModelPart = False

        if OutputMode == "WriteDesignSurface":
            self.WriteDesignSurface = True
        elif OutputMode == "WriteOptimizationModelPart":
            if self.OptimizationModelPart.NumberOfElements() == 0:
                raise NameError("Output of optimization model part in VKT-format requires definition of elements. No elements are given in current mdpa! You may change the design output mode.")
            self.WriteOptimizationModelPart = True
        else:
            raise NameError("The following design output mode is not defined within a VKT output (name may be misspelled): " + OutputMode)

    # --------------------------------------------------------------------------
    def __CreateVTKIO( self ):
        ResultsDirectory = self.OutputSettings["output_directory"].GetString()
        DesignHistoryFilename = self.OutputSettings["design_history_filename"].GetString()
        DesignHistoryFilenameWithPath =  ResultsDirectory+"/"+DesignHistoryFilename

        NodalResults = self.OutputSettings["nodal_results"]

        if self.WriteDesignSurface:
            self.VtkIO = VTKFileIO( self.DesignSurface, DesignHistoryFilenameWithPath, "WriteConditionsOnly", NodalResults )
        elif self.WriteOptimizationModelPart:
            self.VtkIO = VTKFileIO( self.OptimizationModelPart, DesignHistoryFilenameWithPath, "WriteElementsOnly", NodalResults )

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.VtkIO.InitializeLogging()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        self.VtkIO.LogNodalResults( optimizationIteration )

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        pass

# ==============================================================================
