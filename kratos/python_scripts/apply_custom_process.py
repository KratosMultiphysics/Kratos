from KratosMultiphysics import *
import python_process

def Factory(settings, Model):
    """Return an instance of the process.

        Keyword arguments:
        settings -- json parameters containing valid "Parameters" entry
        Model -- dictionary mapping string names to model parts
    """
    params = settings["Parameters"]
    model_part = Model.get(params["model_part_name"].GetString(), "model part not found")

    return ApplyCustomProcess(model_part, params)

class ApplyCustomProcess(python_process.PythonProcess):
    """Execute json list of strings as python code for each process function.

    Public instance variables:
    self.model_part

    Example:
    param = KratosMultiphysics.Parameters('''
    {
      "execute_initialize" : [
        "self.vel_x = 1.0",
        "for node in self.model_part.Nodes:",
        "    node.SetSolutionStepValue(VELOCITY_X, self.vel_x)"
      ],
      "execute_initialize_solution_step" : [
        "self.vel_x += 1.0",
        "for node in self.model_part.Nodes:",
        "    node.SetSolutionStepValue(VELOCITY_X, self.vel_x)"
      ]
    }
    '''
    custom_process = ApplyCustomProcess(model_part, param)
    """
    def __init__(self, model_part, settings):
        """Store model_part as instance variable and convert json lists to multiline strings.

        Each json list of strings forms the body of a process member function.

        Keyword arguments:
        model_part -- a valid model part
        settings -- json parameters
        """
        python_process.PythonProcess.__init__(self)

        default_settings = Parameters("""
            {
                "model_part_name" : "PROVIDE_VALID_MODEL_PART_NAME",
                "execute_initialize" : [],
                "execute_before_solution_loop" : [],
                "execute_initialize_solution_step" : [],
                "execute_finalize_solution_step" : [],
                "execute_before_output_step" : [],
                "execute_after_output_step" : [],
                "execute_finalize" : []
            }
            """
            )
            
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model_part

        self.execute_initialize = ""
        for i in range(settings["execute_initialize"].size()):
            self.execute_initialize += settings["execute_initialize"][i].GetString() + '\n'

        self.execute_before_solution_loop = "" 
        for i in range(settings["execute_before_solution_loop"].size()):
            self.execute_before_solution_loop += settings["execute_before_solution_loop"][i].GetString() + '\n'
            
        self.execute_initialize_solution_step = "" 
        for i in range(settings["execute_initialize_solution_step"].size()):
            self.execute_initialize_solution_step += settings["execute_initialize_solution_step"][i].GetString() + '\n'

        self.execute_finalize_solution_step = "" 
        for i in range(settings["execute_finalize_solution_step"].size()):
            self.execute_finalize_solution_step += settings["execute_finalize_solution_step"][i].GetString() + '\n'

        self.execute_before_output_step = "" 
        for i in range(settings["execute_before_output_step"].size()):
            self.execute_before_output_step += settings["execute_before_output_step"][i].GetString() + '\n'

        self.execute_after_output_step = "" 
        for i in range(settings["execute_after_output_step"].size()):
            self.execute_after_output_step += settings["execute_after_output_step"][i].GetString() + '\n'

        self.execute_finalize = "" 
        for i in range(settings["execute_finalize"].size()):
            self.execute_finalize += settings["execute_finalize"][i].GetString() + '\n'

    def ExecuteInitialize(self):
        exec(self.execute_initialize)
    
    def ExecuteBeforeSolutionLoop(self):
        exec(self.execute_before_solution_loop)
    
    def ExecuteInitializeSolutionStep(self):
        exec(self.execute_initialize_solution_step)
    
    def ExecuteFinalizeSolutionStep(self):
        exec(self.execute_finalize_solution_step)
    
    def ExecuteBeforeOutputStep(self):
        exec(self.execute_before_output_step)
    
    def ExecuteAfterOutputStep(self):
        exec(self.execute_after_output_step)
        
    def ExecuteFinalize(self):
        exec(self.execute_finalize)
