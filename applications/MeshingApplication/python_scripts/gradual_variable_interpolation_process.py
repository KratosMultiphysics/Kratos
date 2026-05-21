import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMA

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return GradualVariableInterpolationProcess(model, settings["Parameters"])

class GradualVariableInterpolationProcess(KratosMultiphysics.Process):
    """This class defines a process for gradually interpolating nodal values from one 
    model part to another, called the origin and destination model parts respectively. 
    The rate of interpolation can be controlled through an increment parameter 'alpha_rampup_increment' 
    or by specifying the total number of steps for the full interpolation 'steps_for_rampup'. 
    The list of variables to be interpolated is also user-definable."""
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "origin_model_part_file_name" : "NameOfMDPAfile",
            "destination_model_part_name" : "ModelPartName",
            "interpolation_variables_list": [],
            "constrain_variables": false,
            "alpha_rampup_increment": 0.0,
            "steps_for_rampup": 0.0
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.settings = settings
        self.constrain_variables = self.settings["constrain_variables"].GetBool()
        self.alpha_rampup_increment = self.settings["alpha_rampup_increment"].GetDouble()
        self.steps_for_rampup = self.settings["steps_for_rampup"].GetInt()

        if not self.alpha_rampup_increment and not self.steps_for_rampup:
            KratosMultiphysics.Logger.PrintWarning('GradualVariableInterpolationProcess', 'Neither alpha_rampup_increment nor steps_for_rampup was specified; defaulting to alpha_rampup_increment = 1.0')
            self.alpha_rampup_increment = 1.0
        elif self.alpha_rampup_increment and self.steps_for_rampup:
            KratosMultiphysics.Logger.PrintWarning('GradualVariableInterpolationProcess', 'Both alpha_rampup_increment and steps_for_rampup were specified; using alpha_rampup_increment and ignoring steps_for_rampup')
        elif self.alpha_rampup_increment:
            if self.alpha_rampup_increment < 0 or self.alpha_rampup_increment > 1:
                raise Exception("'alpha_rampup_increment' must be in the interval (0, 1].")
        elif self.steps_for_rampup:
            if self.steps_for_rampup <= 0:
                raise Exception("'steps_for_rampup' must be a positive integer.")
            else:
                self.alpha_rampup_increment = 1.0 / self.steps_for_rampup

        self.alpha = 0.0
        self.expected_alpha = 1.0
        self.interpolation_variables_list = self.settings["interpolation_variables_list"].GetStringArray()

    def ExecuteInitialize(self):
        origin_model_part_file_name = self.settings["origin_model_part_file_name"].GetString()
        destination_model_part_name = self.settings["destination_model_part_name"].GetString()
        self.destination_model_part = self.model.GetModelPart(destination_model_part_name)
        
        # Import Origin Model Part 
        self.origin_model_part = self.model.CreateModelPart("OriginModelPart")
        for variable in self.destination_model_part.GetHistoricalVariablesNames():
            if KratosMultiphysics.KratosGlobals.HasVariable(variable):
                self.origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))
            else:
                raise Exception(f"Variable '{variable}' does not exist or is not valid in KratosMultiphysics.")
                
        model_part_io = KratosMultiphysics.ModelPartIO(origin_model_part_file_name)
        model_part_io.ReadModelPart(self.origin_model_part)

        self.domain_size = self.destination_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)

        KratosMA.GradualVariableInterpolationUtility.InitializeInterpolationAndConstraints(self.origin_model_part, self.destination_model_part, self.interpolation_variables_list, self.alpha_rampup_increment, self.domain_size, self.constrain_variables)

    def ExecuteInitializeSolutionStep(self):
        self.old_alpha = self.alpha
        if self.old_alpha == 0.0:
            self.old_alpha = self.alpha_rampup_increment
        self.alpha += self.alpha_rampup_increment
        if self.alpha <= self.expected_alpha:
            KratosMA.GradualVariableInterpolationUtility.UpdateSolutionStepVariables(self.destination_model_part, self.interpolation_variables_list, self.alpha, self.old_alpha, self.constrain_variables)