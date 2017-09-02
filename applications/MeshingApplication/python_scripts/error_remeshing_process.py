from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
import os
import math
from json_utilities import *
import json
KratosMultiphysics.CheckForPreviousImport()


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MmgProcess(Model, settings["Parameters"])


class MmgProcess(KratosMultiphysics.Process):

    def __init__(self,Model,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                          : 0,
            "filename"                         : "out",
            "model_part_name"                  : "MainModelPart",
            "strategy"                         : "LevelSet",
            "level_set_strategy_parameters"              :{
                "scalar_variable"                  : "DISTANCE",
                "gradient_variable"                : "DISTANCE_GRADIENT"
            },
            "framework"                            : "Eulerian",
            "internal_variables_parameters"        :
            {
                "allocation_size"                      : 1000, 
                "bucket_size"                          : 4, 
                "search_factor"                        : 2, 
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" :[]
            },
            "hessian_strategy_parameters"              :{
                "metric_variable"                  : ["DISTANCE"],
                "interpolation_error"              : 0.04,
                "mesh_dependent_constant"          : 0.0
            },
            "error_parameters"              :{
                "initial_run"                  : true,
                "interpolation_error"          : 0.04
            },
            "enforce_current"                  : true,
            "initial_step"                     : 1,
            "step_frequency"                   : 0,
            "automatic_remesh"                 : true,
            "automatic_remesh_parameters"      :{
                "automatic_remesh_type"            : "Ratio",
                "min_size_ratio"                   : 1.0,
                "max_size_ratio"                   : 3.0,
                "refer_type"                       : "Mean",
                "min_size_current_percentage"      : 50.0,
                "max_size_current_percentage"      : 98.0
            },
            "initial_remeshing"                : true,
            "fix_contour_model_parts"          : [],
            "minimal_size"                     : 0.1,
            "maximal_size"                     : 10.0,
            "anisotropy_remeshing"             : true,
            "anisotropy_parameters":{
                "hmin_over_hmax_anisotropic_ratio" : 0.01,
                "boundary_layer_max_distance"      : 1.0,
                "boundary_layer_min_size_ratio"    : 2.0,
                "interpolation"                    : "Linear"
            },
            "save_external_files"              : false,
            "max_number_of_searchs"            : 1000,
            "echo_level"                       : 3
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.Model= Model
        self.model_part_name = self.params["model_part_name"].GetString()
        self.dim = self.Model[self.model_part_name].ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.params = params

        self.enforce_current = self.params["enforce_current"].GetBool()

        self.initial_remeshing = self.params["initial_remeshing"].GetBool()

    def ExecuteInitialize(self):
        if (self.dim == 2):
            self.initialize_metric = MeshingApplication.MetricFastInit2D(self.Model[self.model_part_name])
        else:
            self.initialize_metric = MeshingApplication.MetricFastInit3D(self.Model[self.model_part_name])
            
        self.initialize_metric.Execute()

        self._CreateMetricsProcess()

        mmg_parameters = KratosMultiphysics.Parameters("""{}""")
        mmg_parameters.AddValue("filename",self.params["filename"])
        mmg_parameters.AddValue("framework",self.params["framework"])
        mmg_parameters.AddValue("internal_variables_parameters",self.params["internal_variables_parameters"])
        mmg_parameters.AddValue("save_external_files",self.params["save_external_files"])
        mmg_parameters.AddValue("max_number_of_searchs",self.params["max_number_of_searchs"])
        mmg_parameters.AddValue("echo_level",self.params["echo_level"])
        if (self.dim == 2):
            self.MmgProcess = MeshingApplication.MmgProcess2D(self.Model[self.model_part_name], mmg_parameters)
        else:
            self.MmgProcess = MeshingApplication.MmgProcess3D(self.Model[self.model_part_name], mmg_parameters)

        if (self.initial_remeshing == True):
            self._ExecuteRefinement()

    def ExecuteBeforeSolutionLoop(self):
        self.step = 0

    def ExecuteInitializeSolutionStep(self):
        pass

                                
    def ExecuteFinalizeSolutionStep(self):
        if (self.params["error_parameters"]["initial_run"].GetBool() == True):
            self.params["error_parameters"]["initial_run"].SetBool(False)
            self._ExecuteRefinement()
        

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def _CreateMetricsProcess(self):
        self.MetricsProcess = []
        spr_parameters = KratosMultiphysics.Parameters("""{}""")
        spr_parameters.AddValue("minimal_size",self.params["minimal_size"])
        spr_parameters.AddValue("maximal_size",self.params["maximal_size"])
        spr_parameters.AddValue("error",self.params["error_parameters"]["interpolation_error"])
            
        if (self.dim == 2):
            self.MetricsProcess.append(MeshingApplication.ComputeSPRErrorSolMetricProcess2D(
                self.Model[self.model_part_name],
                spr_parameters))
        else:
            self.MetricsProcess.append(MeshingApplication.ComputeSPRErrorSolMetricProcess3D(
                self.Model[self.model_part_name],
                spr_parameters))                        

    def _ExecuteRefinement(self):

        # Initialize metric
        self.initialize_metric.Execute()

        print("Calculating the metrics")
        # Execute metric computation
        for metric_process in self.MetricsProcess:
            metric_process.Execute()

        print("Remeshing")
        self.MmgProcess.Execute()

        # We need to set that the model part has been modified (later on we will act in consequence)
        self.Model[self.model_part_name].Set(KratosMultiphysics.MODIFIED, True)

        print("Remesh finished")