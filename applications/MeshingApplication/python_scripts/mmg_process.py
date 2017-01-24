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
            "input_file_name"                  : "",
            "model_part_name"                  : "MainModelPart",
            "scalar_variable"                  : "DISTANCE",
            "gradient_variable"                : "DISTANCE_GRADIENT",
            "step_frequency"                   : 0,
            "minimal_size"                     : 0.1,
            "hmin_over_hmax_anisotropic_ratio" : 0.01,
            "boundary_layer_max_distance"      : 1.0,
            "interpolation"                    : "Constant",
            "save_external_files"              : false,
            "initialize_nodal_value"           : 1,
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

        self.input_file_name = self.params["input_file_name"].GetString()
        self.minimal_size = self.params["minimal_size"].GetDouble()
        self.hmin_over_hmax_anisotropic_ratio = self.params["hmin_over_hmax_anisotropic_ratio"].GetDouble()
        self.boundary_layer_max_distance = self.params["boundary_layer_max_distance"].GetDouble()
        self.interpolation = self.params["interpolation"].GetString()
        self.step_frequency = self.params["step_frequency"].GetInt()
        self.save_external_files = self.params["save_external_files"].GetBool()
        
        self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["scalar_variable"].GetString() )
        self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["gradient_variable"].GetString() )
        self.initialize_nodal_value = self.params["initialize_nodal_value"].GetInt() 
        self.max_number_of_searchs = self.params["max_number_of_searchs"].GetInt() 
        self.echo_level = self.params["echo_level"].GetInt() 
        
    def ExecuteInitialize(self):          
        self._CreateGradientProcess()
        
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.Model[self.model_part_name])
        
        if (self.dim == 2):
            self.MmgUtility = MeshingApplication.MmgUtility2D(self.input_file_name, self.echo_level)
        else:
            self.MmgUtility = MeshingApplication.MmgUtility3D(self.input_file_name, self.echo_level)
        
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

        print("Preparing the solution and mesh information")
        self.MmgUtility.ComputeExistingModelPart(
            self.Model[self.model_part_name], 
            self.scalar_variable, 
            self.gradient_variable, 
            self.minimal_size, 
            self.hmin_over_hmax_anisotropic_ratio, 
            self.boundary_layer_max_distance, 
            self.interpolation, 
            self.save_external_files
            )
        
        print("Remeshing")
        self.MmgUtility.Execute(self.Model[self.model_part_name], self.save_external_files, self.max_number_of_searchs)
        
        # TODO: Check this
        #if (self.initialize_nodal_value == 1): # Density and kinematic viscosity
            ## Read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
            #for el in self.Model[self.model_part_name].Elements:
                #rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                #kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
                #break

            #KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.Model[self.model_part_name].Nodes)
            #KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.Model[self.model_part_name].Nodes)
        
        self.local_gradient.Execute() # Recalculate gradient after remeshing
