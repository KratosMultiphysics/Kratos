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

# Kratos Core and Apps
import KratosMultiphysics as KM

# Additional imports
from .design_logger_base import DesignLogger
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

# ==============================================================================
class DesignLoggerVTK( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, model_part_controller, optimization_settings ):
        self.output_settings = optimization_settings["output"]
        minimal_vtk_settings = KM.Parameters("""
        {
            "name"       : "vtk",
            "vtk_parameters" : {
                "model_part_name"        : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "output_sub_model_parts" : false,
                "folder_name"            : "Optimization_Results"
            }
        }""")

        output_format = self.output_settings["output_format"]
        if not output_format.Has("vtk_parameters"):
            output_format.AddValue("vtk_parameters", minimal_vtk_settings["vtk_parameters"])
        else:
            if output_format["vtk_parameters"].Has("model_part_name"):
                KM.Logger.PrintWarning("ShapeOpt::DesignLoggerVTK", "vtk output parameter `model_part_name` will be overwritten!")
            if output_format["vtk_parameters"].Has("folder_name"):
                KM.Logger.PrintWarning("ShapeOpt::DesignLoggerVTK", "vtk output parameter `folder_name` will be overwritten!")

        output_format["vtk_parameters"].ValidateAndAssignDefaults(minimal_vtk_settings["vtk_parameters"])

        self.model = model_part_controller.GetModel()
        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.vtk_io = self.__CreateVTKIO()

    # --------------------------------------------------------------------------
    def __CreateVTKIO( self ):

        vtk_parameters = self.output_settings["output_format"]["vtk_parameters"]
        output_mode = self.output_settings["design_output_mode"].GetString()

        nodal_results = self.output_settings["nodal_results"]
        vtk_parameters.AddValue("nodal_solution_step_data_variables", nodal_results)

        if output_mode == "WriteDesignSurface":
            vtk_parameters["model_part_name"].SetString(self.design_surface.FullName())
        elif output_mode == "WriteOptimizationModelPart":
            vtk_parameters["model_part_name"].SetString(self.optimization_model_part.FullName())
        else:
            raise NameError("The following design output mode is not defined within a VTK output (name may be misspelled): " + output_mode)

        vtk_parameters["folder_name"].SetString(self.output_settings["output_directory"].GetString())

        return VtkOutputProcess(self.model, vtk_parameters)

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.vtk_io.ExecuteInitialize()
        self.vtk_io.ExecuteBeforeSolutionLoop()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        OriginalTime = self.optimization_model_part.ProcessInfo[KM.TIME]
        self.optimization_model_part.ProcessInfo[KM.TIME] = optimizationIteration

        self.vtk_io.ExecuteInitializeSolutionStep()
        if(self.vtk_io.IsOutputStep()):
            self.vtk_io.PrintOutput()
        self.vtk_io.ExecuteFinalizeSolutionStep()

        self.optimization_model_part.ProcessInfo[KM.TIME] = OriginalTime

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        self.vtk_io.ExecuteFinalize()

# ==============================================================================
