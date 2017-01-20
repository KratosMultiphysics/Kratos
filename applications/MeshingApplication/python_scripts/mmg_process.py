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
            "mmg_location"            : "~/mmg/bin/mmg3d_O3",
            "input_file_name"         : "",
            "output_file_name"        : "",
            "model_part_name"         : "MainModelPart",
            "scalar_variable"         : "DISTANCE",
            "gradient_variable"       : "DISTANCE_GRADIENT",
            "step_frequency"          : 0,
            "elementary_length"       : 0.1,
            "initial_alpha_parameter" : 0.01,
            "distance_threshold"      : 1.0,
            "interpolation"           : "Constant",
            "save_external_files"     : true,
            "initialize_nodal_value"  : 1,
            "max_number_of_searchs"   : 1000
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
        self.elementary_length = self.params["elementary_length"].GetDouble()
        self.initial_alpha_parameter = self.params["initial_alpha_parameter"].GetDouble()
        self.distance_threshold = self.params["distance_threshold"].GetDouble()
        self.interpolation = self.params["interpolation"].GetString()
        self.mmg_location = self.params["mmg_location"].GetString()
        self.output_file_name = self.params["output_file_name"].GetString()
        self.step_frequency = self.params["step_frequency"].GetInt()
        self.save_external_files = self.params["save_external_files"].GetBool()
        
        self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["scalar_variable"].GetString() )
        self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["gradient_variable"].GetString() )
        self.initialize_nodal_value = self.params["initialize_nodal_value"].GetInt() 
        self.max_number_of_searchs = self.params["max_number_of_searchs"].GetInt() 
        
    def ExecuteInitialize(self):  
        # LEGACY MMG
        ## Basic variables necessary for reading the model_part
        #self.Model[self.model_part_name].AddNodalSolutionStepVariable(self.scalar_variable)
        #self.Model[self.model_part_name].AddNodalSolutionStepVariable(self.gradient_variable)
        #self.Model[self.model_part_name].AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        #self.Model[self.model_part_name].AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VOLUME)
        
        #KratosMultiphysics.ModelPartIO(self.input_file_name).ReadModelPart(self.Model[self.model_part_name])
        
        self._CreateGradientProcess()
        
        if (self.dim == 2):
            self.MmgUtility = MeshingApplication.MmgUtility2D(self.input_file_name)
        else:
            self.MmgUtility = MeshingApplication.MmgUtility3D(self.input_file_name)
        
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
            
        for node in self.Model[self.model_part_name].Nodes:
            node.SetValue(KratosMultiphysics.NODAL_VOLUME, 0.0)
            
        for elem in self.Model[self.model_part_name].Elements:
            area = elem.GetArea()# Area or volume, depending triangle or tetraedra
            if (len(elem.GetNodes()) == 3):
                L = 2.0 * area**(0.5)
            elif (len(elem.GetNodes()) == 4):
                L = 2.0396489026555 * (area)**(1.0/3.0)
            else:
                L = self.elementary_length
                
            for node in elem.GetNodes():
                aux = node.GetValue(KratosMultiphysics.NODAL_VOLUME)
                node.SetValue(KratosMultiphysics.NODAL_VOLUME, aux + L * area )
                
        for node in self.Model[self.model_part_name].Nodes:
            aux = node.GetValue(KratosMultiphysics.NODAL_VOLUME)
            nodal_area = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA, 0)
            if nodal_area > 0.0:
                node.SetValue(KratosMultiphysics.NODAL_VOLUME, aux/nodal_area)

        print("Preparing the solution and mesh information")
        self.MmgUtility.ComputeExistingModelPart(
            self.Model[self.model_part_name], 
            self.scalar_variable, 
            self.gradient_variable, 
            self.elementary_length, 
            self.initial_alpha_parameter, 
            self.distance_threshold, 
            self.interpolation, 
            self.save_external_files
            )
        
        print("Remeshing")
        self.MmgUtility.Execute(self.Model[self.model_part_name], self.save_external_files, self.max_number_of_searchs)
        
        if (self.initialize_nodal_value == 1): # Density and kinematic viscosity
            # Read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
            for el in self.Model[self.model_part_name].Elements:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
                break

            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.Model[self.model_part_name].Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.Model[self.model_part_name].Nodes)
        
        self.local_gradient.Execute() # Recalculate gradient after remeshing
        
        # LEGACY MMG 
        #import write_mmg_mesh as write
        ##write.WriteMmgFile(self.input_file_name, self.elementary_length,self.initial_alpha_parameter, self.distance_threshold, self.model_part_name, self.interpolation)
        #write.WriteMmgFile(self.input_file_name, self.elementary_length,self.initial_alpha_parameter, self.distance_threshold, self.interpolation, self.model_part_name, self.Model[self.model_part_name])
        
        ## Call to the library
        #os.system(self.mmg_location+" -in "+ self.input_file_name +".mesh -out "+ self.output_file_name +".mesh")
        
        #import parse_mmg_file as parse 
        #self.Model[self.model_part_name] = parse.ParseMmgFile(self.output_file_name+".mesh", self.model_part_name)
