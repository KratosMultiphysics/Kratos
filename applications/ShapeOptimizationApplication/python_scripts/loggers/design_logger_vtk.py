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
from KratosMultiphysics.ShapeOptimizationApplication.loggers.design_logger_base import DesignLogger
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
                "file_format"            : "binary",
                "output_sub_model_parts" : false,
                "output_path"            : "Optimization_Results"
            }
        }""")

        output_format = self.output_settings["output_format"]
        if not output_format.Has("vtk_parameters"):
            output_format.AddValue("vtk_parameters", minimal_vtk_settings["vtk_parameters"])
        else:
            if output_format["vtk_parameters"].Has("model_part_name"):
                KM.Logger.PrintWarning("ShapeOpt::DesignLoggerVTK", "vtk output parameter `model_part_name` will be overwritten!")
            if output_format["vtk_parameters"].Has("output_path"):
                KM.Logger.PrintWarning("ShapeOpt::DesignLoggerVTK", "vtk output parameter `output_path` will be overwritten!")

        output_format["vtk_parameters"].ValidateAndAssignDefaults(minimal_vtk_settings["vtk_parameters"])

        self.model = model_part_controller.GetModel()
        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surfaces = model_part_controller.GetDesignSurfaces()

        self.vtk_ios = self.__CreateVTKIO()

    # --------------------------------------------------------------------------
    def __CreateVTKIO( self ):
        vtk_ios = {}
        for name, design_surface in self.design_surfaces.items():
            vtk_parameters = self.output_settings["output_format"]["vtk_parameters"]
            output_mode = self.output_settings["design_output_mode"].GetString()

            nodal_results = self.output_settings["nodal_results"]
            vtk_parameters.AddValue("nodal_solution_step_data_variables", nodal_results)

            if output_mode == "write_design_surface":
                vtk_parameters["model_part_name"].SetString(design_surface.FullName())
            elif output_mode == "write_optimization_model_part":
                vtk_parameters["model_part_name"].SetString(self.optimization_model_part.FullName())
            else:
                raise NameError("The following design output mode is not defined within a VTK output (name may be misspelled): " + output_mode)

            vtk_parameters["output_path"].SetString(self.output_settings["output_directory"].GetString())

            vtk_ios[name] = VtkOutputProcess(self.model, vtk_parameters)

        return vtk_ios

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        for vtk_io in self.vtk_ios.values():
            vtk_io.ExecuteInitialize()
            vtk_io.ExecuteBeforeSolutionLoop()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        OriginalTime = self.optimization_model_part.ProcessInfo[KM.TIME]
        self.optimization_model_part.ProcessInfo[KM.TIME] = optimizationIteration

        for vtk_io in self.vtk_ios.values():
            vtk_io.ExecuteInitializeSolutionStep()
            ## Required for stage wise optimization.
            #if(vtk_io.IsOutputStep()):
            vtk_io.PrintOutput()
            vtk_io.ExecuteFinalizeSolutionStep()

        self.optimization_model_part.ProcessInfo[KM.TIME] = OriginalTime

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        for vtk_io in self.vtk_ios.values():
            vtk_io.ExecuteFinalize()

# ==============================================================================
