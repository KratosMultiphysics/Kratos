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
import KratosMultiphysics

# Import logger base classes
from design_logger_base import DesignLogger

# For VTK output
import KratosMultiphysics.vtk_output_process as vtk_output_process


# ==============================================================================
class DesignLoggerVTK( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, model_part_controller, optimization_settings ):
        self.output_settings = optimization_settings["output"]
        minimal_vtk_settings = KratosMultiphysics.Parameters("""
        {
            "name"              : "vtk",
            "Parameters" : {
                "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "output_sub_model_parts"             : false,
                "folder_name"                        : "Optimization_Results"
            }
        }""")

        output_format = self.output_settings["output_format"]
        if not output_format.Has("Parameters"):
            output_format.AddValue("Parameters", minimal_vtk_settings["Parameters"])
        else:
            if output_format["Parameters"].Has("model_part_name"):
                print("WARNING:: vtk output parameter `model_part_name` will be overwritten!")
            if output_format["Parameters"].Has("output_sub_model_parts"):
                print("WARNING:: vtk output parameter `output_sub_model_parts` will be overwritten!")
            if output_format["Parameters"].Has("folder_name"):
                print("WARNING:: vtk output parameter `folder_name` will be overwritten!")

        output_format["Parameters"].ValidateAndAssignDefaults(minimal_vtk_settings["Parameters"])

        self.model = model_part_controller.GetModel()
        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.vtk_io = self.__CreateVTKIO()

    # --------------------------------------------------------------------------
    def __CreateVTKIO( self ):

        vtk_parameters = self.output_settings["output_format"]["Parameters"]

        nodal_results = self.output_settings["nodal_results"]
        vtk_parameters.AddValue("nodal_solution_step_data_variables", nodal_results)

        output_mode = self.output_settings["design_output_mode"].GetString()
        if output_mode == "WriteDesignSurface":
            vtk_parameters["model_part_name"].SetString(self.design_surface.Name)
        elif output_mode == "WriteOptimizationModelPart":
            vtk_parameters["model_part_name"].SetString(self.optimization_model_part.Name)
        else:
            raise NameError("The following design output mode is not defined within a VTK output (name may be misspelled): " + output_mode)

        vtk_parameters["folder_name"].SetString(self.output_settings["output_directory"].GetString())

        return vtk_output_process.Factory(self.output_settings["output_format"], self.model)

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.vtk_io.ExecuteInitialize()
        self.vtk_io.ExecuteBeforeSolutionLoop()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        OriginalTime = self.optimization_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.optimization_model_part.ProcessInfo[KratosMultiphysics.TIME] = optimizationIteration

        self.vtk_io.ExecuteInitializeSolutionStep()
        if(self.vtk_io.IsOutputStep()):
            self.vtk_io.PrintOutput()
        self.vtk_io.ExecuteFinalizeSolutionStep()

        self.optimization_model_part.ProcessInfo[KratosMultiphysics.TIME] = OriginalTime

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        self.vtk_io.ExecuteFinalize()

# ==============================================================================
