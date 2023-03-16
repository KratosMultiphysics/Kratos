import KratosMultiphysics
import numpy as np
import KratosMultiphysics.MeshingApplication as KratosMA
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignTemperatureProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignTemperatureProcess(KratosMultiphysics.Process):
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
            "computing_model_part_name" : "ComputingModelPartName"
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model = model
        self.settings = settings
        self.temperature = []
        self.velocity_x = []
        self.velocity_y = []

    def ExecuteInitialize(self):
        pass
        

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
        # mp = self.model.GetModelPart("FluidModelPart")
        # for node in mp.Nodes:
        #     if node.Id==2155:
        #         debug = True
        #     node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 1e-2*node.Y)
        #     node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, -1e-2*node.X)
        #     node.SetValue(KratosMultiphysics.VELOCITY_X, 1e-2*node.Y)
        #     node.SetValue(KratosMultiphysics.VELOCITY_Y, -1e-2*node.X)
        #     node.Fix(KratosMultiphysics.VELOCITY_X)
        #     node.Fix(KratosMultiphysics.VELOCITY_Y)

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
        self.computing_model_part_name = self.settings["computing_model_part_name"].GetString()
        mp = self.model.GetModelPart(self.computing_model_part_name)
        # self.computing_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        for node in mp.Nodes:
            if node.Id==2155:
                print(node.X, node.Y, node.GetSolutionStepValue(KratosMultiphysics.VELOCITY), node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
                self.temperature.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
                self.velocity_x.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
                self.velocity_y.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
        debug = True
        # pass

    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        np.save("../velocity_x_co.npy", self.velocity_x)
        np.save("../velocity_y_co.npy", self.velocity_y)
        np.save("../temperature_co.npy", self.temperature)