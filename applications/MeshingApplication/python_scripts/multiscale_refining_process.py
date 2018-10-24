from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MultiscaleRefiningProcess(Model, settings["Parameters"])

class MultiscaleRefiningProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "main_model_part_name"            : "MainModelPart",
            "visualization_model_part_name"   : "VisualizationModelPart",
            "current_subscale"                : 0,
            "maximum_number_of_subscales"     : 4,    
            "echo_level"                      : 0,
            "advanced_configuration"          : {
                "echo_level"                      : 0,
                "number_of_divisions_at_subscale" : 2,
                "subscale_interface_base_name"    : "refined_interface",      
                "subscale_boundary_condition"     : "Condition2D2N"
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.model = Model

        self.echo_level = self.settings['echo_level'].GetInt()
        self.current_subscale = self.settings['current_subscale'].GetInt()
        self.maximum_number_of_subscales = self.settings['maximum_number_of_subscales'].GetInt()

        # Get the coarse model part name
        self.coarse_model_part_name = self.settings['main_model_part_name'].GetString()
        if (self.current_subscale > 0):
            self.coarse_model_part_name += '_' + str(self.current_subscale)
        
        # Get the coarse model part
        self.coarse_model_part = self.model[self.coarse_model_part_name]

        # Get the visualization model part
        if (self.current_subscale == 0):
            self._InitializeVisualizationModelPart()
        else:
            self.visualization_model_part = self.model[self.settings['visualization_model_part_name'].GetString()]

        # Initialize the refined model part
        if (self.current_subscale < self.maximum_number_of_subscales):
            self._InitializeRefinedModelPart()

    def ExecuteInitialize(self):
        # Create the new subscale process
        self.subscales_utility = MeshingApplication.MultiscaleRefiningProcess(
            self.coarse_model_part,
            self.refined_model_part,
            self.visualization_model_part,
            self.settings["advanced_configuration"])

        if self.echo_level > 0:
            print('The multiscale process is initialized')

        if self.echo_level > 1:
            print(self.model[self.coarse_model_part_name])
            print(self.model[self.refined_model_part_name])

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        self._ExecuteRefinement()

    def ExecuteFinalizeSolutionStep(self):
        self._ExecuteCoarsening()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def Clear(self):
        pass

    def GetCoarseModelPart(self):
        return self.coarse_model_part

    def GetRefinedModelPart(self):
        return self.refined_model_part

    def GetVisualizationModelPart(self):
        return self.visualization_model_part

    def _InitializeVisualizationModelPart(self):
        visualization_model_part_name = self.settings['visualization_model_part_name'].GetString()
        self.visualization_model_part = self.model.CreateModelPart(visualization_model_part_name)
        buffer_size = self.coarse_model_part.GetBufferSize()
        self.visualization_model_part.SetBufferSize(buffer_size)
        MeshingApplication.MultiscaleRefiningProcess.InitializeNewModelPart(self.coarse_model_part, self.visualization_model_part)

    def _InitializeRefinedModelPart(self):
        self.refined_model_part_name = self.settings['main_model_part_name'].GetString() + '_' + str(self.current_subscale + 1)
        self.refined_model_part = self.model.CreateModelPart(self.refined_model_part_name)
        buffer_size = self.coarse_model_part.GetBufferSize()
        self.refined_model_part.SetBufferSize(buffer_size)
        MeshingApplication.MultiscaleRefiningProcess.InitializeNewModelPart(self.coarse_model_part, self.refined_model_part)

    def _ExecuteRefinement(self):
        self.subscales_utility.ExecuteRefinement()

    def _ExecuteCoarsening(self):
        self.subscales_utility.ExecuteCoarsening()