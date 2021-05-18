# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
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
            "refining_condition_parameters"   : {
                "kratos_module"                   : "KratosMultiphysics.ShallowWaterApplication",
                "python_module"                   : "residual_based_refining_condition_process",
                "Parameters"                      : {
                    "model_part_name"                 : "model_part",
                    "error_variable"                  : "RESIDUAL_NORM",
                    "variable_threshold"              : 1e-3,
                    "increase_threshold"              : true,
                    "only_refine_wet_domain"          : true
                }
            },
            "variables_to_apply_fixity"       : [],
            "variables_to_set_at_interface"   : [],
            "variables_to_update_at_coarse"   : [],
            "advanced_configuration"          : {
                "echo_level"                      : 0,
                "number_of_divisions_at_subscale" : 2,
                "subscale_interface_base_name"    : "refined_interface",
                "subscale_boundary_condition"     : "LineCondition2D2N"
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
        self.number_of_divisions_at_subscale = self.settings['advanced_configuration']['number_of_divisions_at_subscale'].GetInt()
        self.number_of_substeps = 2**self.settings['advanced_configuration']['number_of_divisions_at_subscale'].GetInt()

        # Get the coarse and refined model part names
        self.coarse_model_part_name = self.settings['main_model_part_name'].GetString()
        self.refined_model_part_name = self.coarse_model_part_name
        if (self.current_subscale == 0):
            self.coarse_model_part_name += '_0'
            self.refined_model_part_name += '_0'
        else:
            self.coarse_model_part_name += '_' + str(self.current_subscale - 1)
            self.refined_model_part_name += '_' + str(self.current_subscale)

        # Get the model parts
        if (self.current_subscale == 0):
            self.coarse_model_part = self.model.CreateModelPart(self.coarse_model_part_name)
            self.visualization_model_part = self.model.CreateModelPart(self.settings['visualization_model_part_name'].GetString())
            self.refined_model_part = self.coarse_model_part
        else:
            self.coarse_model_part = self.model[self.coarse_model_part_name]
            self.visualization_model_part = self.model[self.settings['visualization_model_part_name'].GetString()]
            self.refined_model_part = self.model.CreateModelPart(self.refined_model_part_name)
            self.refined_model_part.ProcessInfo[MeshingApplication.SUBSCALE_INDEX] = self.current_subscale

        self.variables_to_apply_fixity = self._GenerateVariableListFromInput(self.settings["variables_to_apply_fixity"])
        self.variables_to_set_at_interface = self._GenerateVariableListFromInput(self.settings["variables_to_set_at_interface"])
        self.variables_to_update_at_coarse = self._GenerateVariableListFromInput(self.settings["variables_to_update_at_coarse"])

        kratos_module_name = self.settings["refining_condition_parameters"]["kratos_module"].GetString()
        python_module_name = self.settings["refining_condition_parameters"]["python_module"].GetString()
        full_module_name = kratos_module_name + "." + python_module_name
        python_module = __import__(full_module_name, fromlist=[python_module_name])
        self.settings["refining_condition_parameters"]["Parameters"]["model_part_name"].SetString(self.coarse_model_part_name)
        self.refining_condition_process = python_module.Factory(self.settings["refining_condition_parameters"], self.model)

    def PrepareModelPart(self):
        # Initialize the corresponding model part
        if self.current_subscale == 0:
            self._InitializeVisualizationModelPart()
        else:
            self._InitializeRefinedModelPart()
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

    def ExecuteInitialize(self):
        # This is equivalent to the ExecuteBeforeSolutionLoop
        self.ExecuteBeforeSolutionLoop()

    def ExecuteBeforeSolutionLoop(self):
        if self.current_subscale > 0:
            self._EvaluateCondition()
            self._ExecuteRefinement()
            self._ExecuteCoarsening()
            self._ApplyFixityAtInterface(True)

    def ExecuteInitializeSolutionStep(self):
        if self.current_subscale == 0:
            time = self.coarse_model_part.ProcessInfo[KratosMultiphysics.TIME]
            step = self.coarse_model_part.ProcessInfo[KratosMultiphysics.STEP]
            self.visualization_model_part.ProcessInfo[KratosMultiphysics.TIME] = time
            self.visualization_model_part.ProcessInfo[KratosMultiphysics.STEP] = step
        else:
            self._TransferSubstepToRefinedInterface()

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        if self.current_subscale > 0:
            self._ApplyFixityAtInterface(False)
            self._UpdateVariablesAtCoarseModelPart()

    def Clear(self):
        pass

    def GetCoarseModelPart(self):
        return self.coarse_model_part

    def GetRefinedModelPart(self):
        return self.refined_model_part

    def GetVisualizationModelPart(self):
        return self.visualization_model_part

    def _InitializeVisualizationModelPart(self):
        buffer_size = self.coarse_model_part.GetBufferSize()
        self.visualization_model_part.SetBufferSize(buffer_size)
        MeshingApplication.MultiscaleRefiningProcess.InitializeVisualizationModelPart(self.coarse_model_part, self.visualization_model_part)

    def _InitializeRefinedModelPart(self):
        buffer_size = self.coarse_model_part.GetBufferSize()
        self.refined_model_part.SetBufferSize(buffer_size)
        MeshingApplication.MultiscaleRefiningProcess.InitializeRefinedModelPart(self.coarse_model_part, self.refined_model_part)

    def _ExecuteRefinement(self):
        self.subscales_utility.ExecuteRefinement()

    def _ExecuteCoarsening(self):
        self.subscales_utility.ExecuteCoarsening()

    def _EvaluateCondition(self):
        self.refining_condition_process.Execute()

    def _ApplyFixityAtInterface(self, state):
        for variable in self.variables_to_apply_fixity:
            self.subscales_utility.FixRefinedInterface(variable, state)

    def _TransferSubstepToRefinedInterface(self):
        substep_fraction = self.refined_model_part.ProcessInfo[KratosMultiphysics.STEP] / self.number_of_substeps
        for variable in self.variables_to_set_at_interface:
            self.subscales_utility.TransferSubstepToRefinedInterface(variable, substep_fraction)

    def _UpdateVariablesAtCoarseModelPart(self):
        for variable in self.variables_to_update_at_coarse:
            self.subscales_utility.TransferLastStepToCoarseModelPart(variable)

    def _GenerateVariableListFromInput(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]
