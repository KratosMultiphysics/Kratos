import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMA

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return GradualVariableInterpolationProcess(model, settings["Parameters"])

class GradualVariableInterpolationProcess(KratosMultiphysics.Process):
    """This class defines a process for gradually interpolating nodal values from one 
    model part to another, called the origin and destination model parts respectively. 
    The rate of interpolation can be controlled through an increment parameter 'alpha_increment' 
    or by specifying the total number of steps for the full interpolation 'steps_for_full_interpolation'. 
    The list of variables to be interpolated is also user-definable."""
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "origin_model_part_file_name" : "NameOfMDPAfile",
            "destination_model_part_name" : "ModelPartName",
            "interpolation_variables_list": [],
            "constrain_varibles": false,
            "alpha_increment": 0.0,
            "steps_for_full_interpolation": 0.0
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.settings = settings
        self.constrain_variables = self.settings["constrain_varibles"].GetBool()
        self.alpha_increment = self.settings["alpha_increment"].GetDouble()
        self.steps_for_full_interpolation = self.settings["steps_for_full_interpolation"].GetInt()

        if not self.alpha_increment and not self.steps_for_full_interpolation:
            KratosMultiphysics.Logger.PrintWarning('GradualVariableInterpolationProcess', 'Neither alpha_increment nor steps_for_full_interpolation was specified; defaulting to alpha_increment = 1.0')
            self.alpha_increment = 1.0
        elif self.alpha_increment and self.steps_for_full_interpolation:
            KratosMultiphysics.Logger.PrintWarning('GradualVariableInterpolationProcess', 'Both alpha_increment and steps_for_full_interpolation were specified; using alpha_increment and ignoring steps_for_full_interpolation')
        elif self.alpha_increment:
            if self.alpha_increment < 0 or self.alpha_increment > 1:
                raise Exception("'alpha_increment' must be in the interval (0, 1].")
        elif self.steps_for_full_interpolation:
            if self.steps_for_full_interpolation <= 0:
                raise Exception("'steps_for_full_interpolation' must be a positive integer.")
            else:
                self.alpha_increment = 1.0 / self.steps_for_full_interpolation

        self.alpha = 0.0
        self.expected_alpha = 1.0
        self.interpolation_variables_list = [KratosMultiphysics.KratosGlobals.GetVariable(i) for i in self.settings["interpolation_variables_list"].GetStringArray()]

    def ExecuteInitialize(self):
        self.origin_model_part_file_name = self.settings["origin_model_part_file_name"].GetString()
        destination_model_part_name = self.settings["destination_model_part_name"].GetString()
        self.destination_model_part = self.model.GetModelPart(destination_model_part_name)

        self.InterpolateVariables()

    def InterpolateVariables(self):
        #Import Origin Model Part 
        self.origin_model_part = self.model.CreateModelPart("OriginModelPart")
        for variable in self.destination_model_part.GetHistoricalVariablesNames():
            if KratosMultiphysics.KratosGlobals.HasVariable(variable):
                self.origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))

        model_part_io = KratosMultiphysics.ModelPartIO(self.origin_model_part_file_name)
        model_part_io.ReadModelPart(self.origin_model_part) 

        #Set velocity field to origin model part
        for node in self.origin_model_part.Nodes:
            for variable in self.interpolation_variables_list:
                node.SetValue(variable, self.alpha_increment * node.GetSolutionStepValue(variable))

        #Interpolate velocity to destination model part 
        interpolation = KratosMA.NodalValuesInterpolationProcess2D(self.origin_model_part, self.destination_model_part)
        interpolation.Execute()
        if self.constrain_variables:
            for node in self.destination_model_part.Nodes:
                for variable in self.interpolation_variables_list:
                    node.SetSolutionStepValue(variable, node.GetValue(variable))
                    node.Fix(variable)
        else:
            for node in self.destination_model_part.Nodes:
                for variable in self.interpolation_variables_list:
                    node.SetSolutionStepValue(variable, node.GetValue(variable))


    def ExecuteInitializeSolutionStep(self):
        self.old_alpha = self.alpha
        if self.old_alpha == 0.0:
            self.old_alpha = self.alpha_increment
        self.alpha += self.alpha_increment
        if self.alpha <= self.expected_alpha:
            if self.constrain_variables:
                for node in self.destination_model_part.Nodes:
                    for variable in self.interpolation_variables_list:
                        node.SetSolutionStepValue(variable, self.alpha * node.GetSolutionStepValue(variable) / self.old_alpha)
                        node.Fix(variable)
            else:
                for node in self.destination_model_part.Nodes:
                    for variable in self.interpolation_variables_list:
                        node.SetSolutionStepValue(variable, self.alpha * node.GetSolutionStepValue(variable) / self.old_alpha)

