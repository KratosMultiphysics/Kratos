from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
import os
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
            "output_file_name"                  : "",
            "model_part_name"                  : "MainModelPart",
            "strategy"                         : "Level Set",
            "strategy_parameters"              :{
                "scalar_variable"                  : "DISTANCE",
                "gradient_variable"                : "DISTANCE_GRADIENT"
            },
            "step_frequency"                   : 0,
            "minimal_size"                     : 0.1,
            "anisotropy_remeshing"             : true,
            "anisotropy_parameters":{
                "hmin_over_hmax_anisotropic_ratio" : 0.01,
                "boundary_layer_max_distance"      : 1.0,
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
        
        self.strategy = self.params["strategy"].GetString()
        if (self.strategy == "Level Set"):
            self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["strategy_parameters"]["scalar_variable"].GetString() )
            self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["strategy_parameters"]["gradient_variable"].GetString() )
        elif (self.strategy == "Hessian"):
            self.metric_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["metric_variable"].GetString() )
            
        self.minimal_size = self.params["minimal_size"].GetDouble()
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
        self._CreateGradientProcess()
        
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.Model[self.model_part_name])
        
        if (self.dim == 2):
            self.MmgUtility = MeshingApplication.MmgUtility2D(self.output_file_name, self.echo_level)
        else:
            self.MmgUtility = MeshingApplication.MmgUtility3D(self.output_file_name, self.echo_level)
        
        self._ExecuteRefinement()
        
    def ExecuteBeforeSolutionLoop(self):
        self.step = 0
        
    def ExecuteInitializeSolutionStep(self):
        self.step += 1
        if self.step_frequency > 0: # NOTE: Requires to recalculate the scalar value in each node (we interpolate)
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
            
        self.local_gradient.Execute()
        
        self.find_nodal_h.Execute()

        print("Preparing the solution  and mesh information")
        self.MmgUtility.InitializeMeshData(self.Model[self.model_part_name])
        if (self.strategy == "Level Set"):
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
        
        self.local_gradient.Execute() # Recalculate gradient after remeshing
