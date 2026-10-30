import KratosMultiphysics as KM
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA


def Factory(settings, Model):
    """Factory function to create the process.

    Keyword arguments:
    settings -- Kratos parameters containing solver settings.
    Model -- the container of the different model parts.
    """
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    
    # Extract the model part name and parameters
    model_part_name = settings["Parameters"]["model_part_name"].GetString()
    model_part = Model[model_part_name]
    process_info = model_part.ProcessInfo
    
    return CSMA.PrepareLinearConstraintsForQuadraticContactProcess(model_part, process_info)


class PrepareLinearConstraintsForQuadraticContactProcess(KM.Process):
    """This process is used for preparing linear constraints for quadratic contact.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """The default constructor of the class.

        Keyword arguments:
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KM.Process.__init__(self)
        
        # Default parameters
        default_parameters = KM.Parameters("""
        {
            "model_part_name"       : "Structure",
            "help"                  : "This process is used for preparing linear constraints for quadratic contact"
        }
        """)
        
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)
        
        # Get the model part
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.main_model_part = Model[self.model_part_name]
        self.process_info = self.main_model_part.ProcessInfo
        
        # Create the C++ process
        self.cpp_process = CSMA.PrepareLinearConstraintsForQuadraticContactProcess(
            self.main_model_part, 
            self.process_info
        )

    def ExecuteProcess(self):
        """This method calls the C++ ExecuteProcess method."""
        self.cpp_process.ExecuteProcess()

    def ExecuteInitialize(self):
        """This method is executed at the beginning to initialize the process."""
        self.cpp_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """This method is executed before starting the time loop."""
        self.cpp_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """This method is executed in order to initialize the current step."""
        self.cpp_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """This method is executed in order to finalize the current step."""
        self.cpp_process.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        """This method is executed right before the output process computation."""
        self.cpp_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """This method is executed right after the output process computation."""
        self.cpp_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """This method is executed in order to finalize the current computation."""
        self.cpp_process.ExecuteFinalize()

    def Check(self):
        """This function is designed for being called after ExecuteInitialize ONCE
        to verify that the input is correct.
        """
        return self.cpp_process.Check()
