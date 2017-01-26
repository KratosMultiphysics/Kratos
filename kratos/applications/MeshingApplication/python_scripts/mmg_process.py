from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
import os
KratosMultiphysics.CheckForPreviousImport()

import math
def normpdf(x, mean, sd):
    var = float(sd)**2
    pi = 3.1415926
    denom = (2*pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def normvalf(prob, mean, sd):
    var = float(sd)**2
    pi = 3.1415926
    factor = (2*pi*var)**.5
    x = (- math.log(prob * factor)*(2*var))**(0.5) + float(mean)
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
                "metric_variable"                  : "DISPLACEMENT"
            },
            "step_frequency"                   : 0,
            "automatic_remesh"                 : true,
            "automatic_remesh_parameters"      :{
                "automatic_remesh_type"            : "Ratio",
                "refer_type"                       : "Mean",
                "size_ratio"                       : 1.0,
                "size_current_percentage"          : 50
            },
            "minimal_size"                     : 0.1,
            "anisotropy_remeshing"             : true,
            "anisotropy_parameters":{
                "hmin_over_hmax_anisotropic_ratio" : 0.01,
                "boundary_layer_max_distance"      : 1.0,
                "boundary_layer_size_ratio"        : 1.0,
                "interpolation"                    : "Linear"
            },
            "save_external_files"              : false,
            "max_number_of_searchs"            : 1000,
            "echo_level"                       : 3
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        
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
        
        # Calculate NODAL_H
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.Model[self.model_part_name])
        self.find_nodal_h.Execute()
        
        # Calculate the parameters of automatic remeshing
        if (self.params["automatic_remesh"].GetBool() == True):
            import statistics as stat
            nodal_h_values = []
            for node in self.Model[self.model_part_name].Nodes:
                nodal_h_values.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_H))
            
            # NOTE: For mode: https://docs.python.org/3/library/statistics.html
            if (self.params["automatic_remesh_parameters"]["refer_type"].GetString() == "Mean"):
                ref = stat.mean(nodal_h_values)
            elif (self.params["automatic_remesh_parameters"]["refer_type"].GetString() == "Median"):
                ref = stat.median(nodal_h_values)
            
            # Calculate the minimum size
            if (self.params["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Ratio"):
                self.minimal_size = ref * self.params["automatic_remesh_parameters"]["size_ratio"].GetDouble()
            elif (self.params["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Percentage"):
                mean = stat.mean(nodal_h_values)
                stdev = stat.stdev(nodal_h_values)
                prob = self.params["automatic_remesh_parameters"]["size_current_percentage"].GetDouble()
                self.minimal_size = normvalf(prob, mean, stdev)
        else:
            # Manually defined
            self.minimal_size = self.params["minimal_size"].GetDouble()
            
        # Anisotropic remeshing parameters
        self.anisotropy_remeshing = self.params["anisotropy_remeshing"].GetBool()
        if (self.anisotropy_remeshing == True):
            self.hmin_over_hmax_anisotropic_ratio = self.params["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble()
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
            
        if (self.strategy == "LevelSet"):
            # Calculate the gradient
            self.local_gradient.Execute()
        
        # Recalculate NODAL_H
        self.find_nodal_h.Execute()

        print("Preparing the solution  and mesh information")
        self.MmgUtility.InitializeMeshData(self.Model[self.model_part_name])
        if (self.strategy == "LevelSet"):
            if (self.anisotropy_remeshing == True):
                self.MmgUtility.InitializeLevelSetSolData(
                    self.Model[self.model_part_name], 
                    self.minimal_size, 
                    self.scalar_variable, 
                    self.gradient_variable, 
                    self.hmin_over_hmax_anisotropic_ratio, 
                    self.boundary_layer_max_distance, 
                    self.interpolation
                    )
            else:
                self.MmgUtility.InitializeLevelSetSolData( 
                    self.Model[self.model_part_name], 
                    self.minimal_size
                    )
        #elif (self.strategy == "Hessian"):
            #if (self.anisotropy_remeshing == True):
                #self.MmgUtility.InitializeHessianSolData(
                    #self.Model[self.model_part_name], 
                    #self.minimal_size
                    #)
            #else:
                #self.MmgUtility.InitializeHessianSolData( 
                    #self.Model[self.model_part_name], 
                    #self.minimal_size
                    #)
        
        print("Remeshing")
        self.MmgUtility.RemeshModelPart(self.Model[self.model_part_name], self.save_external_files, self.max_number_of_searchs)
        
        if (self.strategy == "LevelSet"):
            self.local_gradient.Execute() # Recalculate gradient after remeshing
