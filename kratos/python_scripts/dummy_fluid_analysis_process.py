import KratosMultiphysics
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DummyFluidAnalysisProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class DummyFluidAnalysisProcess(KratosMultiphysics.Process):
    """This class is a dummy-process that shows how the functions that can be implemented
    in order to customize the behavior

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KratosMultiphysics.Process.__init__(self) # calling the baseclass constructor

        default_settings = KratosMultiphysics.Parameters("""{
            "domain_model_part_name" : "ComputingModelPartName",
            "max_velocity" : 0.0
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model = model
        self.settings = settings

    def ExecuteInitialize(self):
        self.domain_model_part_name = self.settings["domain_model_part_name"].GetString()
        self.domain_model_part = self.model.GetModelPart(self.domain_model_part_name)
        self.max_velocity = self.settings["max_velocity"].GetDouble()
        self.velocity_factor_flag = True
        
        

    def Check(self):
        """ This method verifies that the input is correct

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed just before the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.velocity_factor_flag:
            self.velocity_factor = (np.sin(self.domain_model_part.ProcessInfo[KratosMultiphysics.TIME])**2)
            if self.velocity_factor>0.99:
                self.velocity_factor = 1
                self.velocity_factor_flag = False
        for node in self.domain_model_part.Nodes:
            value = self.max_velocity*self.velocity_factor*10*(-1000*node.Z**2-1000*node.Y**2+0.1)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, value)
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)

    def ExecuteBeforeOutputStep(self):
        """ This method is executed before writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed after writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass