from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
import os
from json_utilities import *
import json
KratosMultiphysics.CheckForPreviousImport()

def linear_interpolation(x, x_list, y_list):
    ind_inf = 0
    ind_sup = -1
    x_inf = x_list[ind_inf]
    x_sup = x_list[ind_sup]
    for i in range(len(x_list)):
        if x_list[i] <= x:
            ind_inf = i
            x_inf = x_list[ind_inf]
        if x_list[-(1 + i)] >= x:
            ind_sup = -(1 + i)
            x_sup = x_list[ind_sup]
    
    if (x_sup - x_inf == 0):
        y = y_list[ind_inf]
    else:
        prop_sup = (x - x_inf)/(x_sup - x_inf)
        prop_inf = 1.0 - prop_sup
        y = y_list[ind_inf] * prop_inf + y_list[ind_sup] * prop_sup 
    
    return y

def normpdf(x, mean, sd):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data =  read_external_json(dir_path+"/normal_distribution.json")
    z = (x-mean)/sd
    z_list    = data["Z"]
    prob_list = data["Prob"]
    if (z > 0):
        prob = linear_interpolation(z, z_list, prob_list)
    else:
        prob = 1.0 - linear_interpolation(-z, z_list, prob_list)
    return prob

def normvalf(prob, mean, sd):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data =  read_external_json(dir_path+"/normal_distribution.json")
    z_list    = data["Z"]
    prob_list = data["Prob"]
    if (prob >= 0.5):
        z = linear_interpolation(prob, prob_list, z_list)
    else:
        z = - linear_interpolation(1.0 - prob, prob_list, z_list)
    x = z * sd + mean
    return x

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MmgProcess(Model, settings["Parameters"])

class MmgProcess(KratosMultiphysics.Process):
  
    def __init__(self,Model,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "output_file_name"                  : "",
            "model_part_name"                  : "MainModelPart",
            "strategy"                         : "LevelSet",
            "level_set_strategy_parameters"              :{
                "scalar_variable"                  : "DISTANCE",
                "gradient_variable"                : "DISTANCE_GRADIENT"
            },
            "hessian_strategy_parameters"              :{
                "metric_variable"                  : "DISPLACEMENT",
                "interpolation_error"              : 1.0e-6,
                "mesh_dependent_constant"          : 0.28125
            },
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
        
        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.RecursivelyValidateAndAssignDefaults(default_parameters)
        
        self.Model= Model
        self.model_part_name = self.params["model_part_name"].GetString()
        self.dim = self.Model[self.model_part_name].ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.params = params

        self.output_file_name = self.params["output_file_name"].GetString()
        
        # Select the remeshing strategy
        self.strategy = self.params["strategy"].GetString()
        if (self.strategy == "LevelSet"):
            self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["level_set_strategy_parameters"]["scalar_variable"].GetString() )
            self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["level_set_strategy_parameters"]["gradient_variable"].GetString() )
        elif (self.strategy == "Hessian"):
            self.metric_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["hessian_strategy_parameters"]["metric_variable"].GetString() )
            self.interpolation_error = self.params["hessian_strategy_parameters"]["interpolation_error"].GetDouble()  
            self.mesh_dependent_constant = self.params["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble()  
        
        # Calculate NODAL_H
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.Model[self.model_part_name])
        self.find_nodal_h.Execute()
        
        # Calculate the parameters of automatic remeshing
        if (self.params["automatic_remesh"].GetBool() == True):
            import statistics as stat
            nodal_h_values = []
            for node in self.Model[self.model_part_name].Nodes:
                nodal_h_values.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_H))
            
            # Calculate the minimum size
            if (self.params["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Ratio"):
                # NOTE: For mode: https://docs.python.org/3/library/statistics.html
                if (self.params["automatic_remesh_parameters"]["refer_type"].GetString() == "Mean"):
                    ref = stat.mean(nodal_h_values)
                elif (self.params["automatic_remesh_parameters"]["refer_type"].GetString() == "Median"):
                    ref = stat.median(nodal_h_values)
                
                self.minimal_size = ref * (self.params["automatic_remesh_parameters"]["min_size_ratio"].GetDouble())
                self.maximal_size = ref * (self.params["automatic_remesh_parameters"]["max_size_ratio"].GetDouble())
            elif (self.params["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Percentage"):
                mean = stat.mean(nodal_h_values)
                stdev = stat.stdev(nodal_h_values)
                prob = (self.params["automatic_remesh_parameters"]["min_size_current_percentage"].GetDouble())/100
                self.minimal_size = normvalf(prob, mean, stdev)
                
                prob = (self.params["automatic_remesh_parameters"]["max_size_current_percentage"].GetDouble())/100
                self.maximal_size = normvalf(prob, mean, stdev)
        else:
            # Manually defined
            self.minimal_size = self.params["minimal_size"].GetDouble()
            self.maximal_size = self.params["maximal_size"].GetDouble()
            
        # Anisotropic remeshing parameters
        self.anisotropy_remeshing = self.params["anisotropy_remeshing"].GetBool()
        if (self.anisotropy_remeshing == True):
            self.hmin_over_hmax_anisotropic_ratio = self.params["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble()
            if (self.params["automatic_remesh"].GetBool() == True):
                self.boundary_layer_max_distance = self.minimal_size * self.params["anisotropy_parameters"]["boundary_layer_min_size_ratio"].GetDouble()
            else:
                self.boundary_layer_max_distance = self.params["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble()
            self.interpolation = self.params["anisotropy_parameters"]["interpolation"].GetString()
            
        self.step_frequency = self.params["step_frequency"].GetInt()
        self.save_external_files = self.params["save_external_files"].GetBool()
        self.max_number_of_searchs = self.params["max_number_of_searchs"].GetInt() 
        self.echo_level = self.params["echo_level"].GetInt() 
        
    def ExecuteInitialize(self): 
        if (self.strategy == "LevelSet"):
            self._CreateGradientProcess()
        
        if (self.dim == 2):
            self.MmgUtility = MeshingApplication.MmgUtility2D(self.output_file_name, self.echo_level)
        else:
            self.MmgUtility = MeshingApplication.MmgUtility3D(self.output_file_name, self.echo_level)
        
        self._ExecuteRefinement()
        
    def ExecuteBeforeSolutionLoop(self):
        self.step = 0
        
    def ExecuteInitializeSolutionStep(self):
        self.step += 1
        if self.step_frequency > 0:
            if self.step == self.step_frequency:
                self._ExecuteRefinement()
                self.step = 0 # Reset
            
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass
    
    def _CreateGradientProcess(self):
        # We compute the scalar value gradient
        if (self.dim == 2):
            self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.Model[self.model_part_name], self.scalar_variable, self.gradient_variable, KratosMultiphysics.NODAL_AREA)
        else: 
            self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.Model[self.model_part_name], self.scalar_variable, self.gradient_variable, KratosMultiphysics.NODAL_AREA)
            
    def _ExecuteRefinement(self):
        
        self._InitializeMetric()
            
        if (self.strategy == "LevelSet"):
            # Calculate the gradient
            self.local_gradient.Execute()
        
        # Recalculate NODAL_H
        self.find_nodal_h.Execute()

        print("Calculating the metrics")
        if (self.dim == 2):
            if (self.anisotropy_remeshing == True):
                MetricsUtility = MeshingApplication.MetricsUtility2D(
                    self.minimal_size,
                    self.hmin_over_hmax_anisotropic_ratio, 
                    self.boundary_layer_max_distance, 
                    self.interpolation
                    )
            else:
                MetricsUtility = MeshingApplication.MetricsUtility2D(self.minimal_size)
        else:
            if (self.anisotropy_remeshing == True):
                MetricsUtility = MeshingApplication.MetricsUtility3D(
                    self.minimal_size,
                    self.hmin_over_hmax_anisotropic_ratio, 
                    self.boundary_layer_max_distance, 
                    self.interpolation
                    )
            else:
                MetricsUtility = MeshingApplication.MetricsUtility3D(self.minimal_size)
            
        if (self.strategy == "LevelSet"):
            MetricsUtility.ComputeLevelSetSolMetric(                    self.Model[self.model_part_name], self.gradient_variable)
  
        elif (self.strategy == "Hessian"):
            for node in self.Model[self.model_part_name].Nodes:
                val = node.GetSolutionStepValue(self.metric_variable, 0)
                break
            if isinstance(val,float): 
                MetricsUtility.ComputeHessianMetric(
                        self.Model[self.model_part_name], 
                        self.maximal_size,
                        self.metric_variable,
                        self.interpolation_error,
                        self.mesh_dependent_constant
                        )
            else:
                components = [KratosMultiphysics.KratosGlobals.GetVariable( self.params["hessian_strategy_parameters"]["metric_variable"].GetString() + "_X" ),KratosMultiphysics.KratosGlobals.GetVariable( self.params["hessian_strategy_parameters"]["metric_variable"].GetString() + "_Y" )]
                if (self.dim == 3):
                    components.append(KratosMultiphysics.KratosGlobals.GetVariable( self.params["hessian_strategy_parameters"]["metric_variable"].GetString() + "_Z" ))
                for comp in components:
                    MetricsUtility.ComputeHessianMetric(
                            self.Model[self.model_part_name], 
                            self.maximal_size,
                            comp,
                            self.interpolation_error,
                            self.mesh_dependent_constant
                            )
        
        print("Remeshing")
        self.MmgUtility.RemeshModelPart(self.Model[self.model_part_name], self.save_external_files, self.max_number_of_searchs)
        
        if (self.strategy == "LevelSet"):
            self.local_gradient.Execute() # Recalculate gradient after remeshing
            
    def _InitializeMetric(self):
        # Initialize metric
        if (self.dim == 2):
            ZeroVector = KratosMultiphysics.Vector(3) 
            ZeroVector[0] = 0.0
            ZeroVector[1] = 0.0
            ZeroVector[2] = 0.0
        else:
            ZeroVector = KratosMultiphysics.Vector(6) 
            ZeroVector[0] = 0.0
            ZeroVector[1] = 0.0
            ZeroVector[2] = 0.0
            ZeroVector[3] = 0.0
            ZeroVector[4] = 0.0
            ZeroVector[5] = 0.0
        for node in self.Model[self.model_part_name].Nodes:
            node.SetValue(MeshingApplication.MMG_METRIC, ZeroVector)
