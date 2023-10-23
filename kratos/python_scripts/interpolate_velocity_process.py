import KratosMultiphysics
import numpy as np
import KratosMultiphysics.MeshingApplication as KratosMA

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InterpolateVelocityProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class InterpolateVelocityProcess(KratosMultiphysics.Process):
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
            "origin_model_part_file_name" : "NameOfMDPAfile",
            "destination_model_part_name" : "ModelPartName",
            "outlet_model_part_name": "ModelPartName",
            "ambient_temperature": 0.0
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model = model
        self.settings = settings

    def ExecuteInitialize(self):
        self.origin_model_part_file_name = self.settings["origin_model_part_file_name"].GetString()
        destination_model_part_name = self.settings["destination_model_part_name"].GetString()
        self.outlet_model_part_name = self.settings["outlet_model_part_name"].GetString()
        self.ambient_temperature = self.settings["ambient_temperature"].GetDouble()
        self.destination_model_part = self.model.GetModelPart(destination_model_part_name)
        

        self.InterpolateVelocityWithMA()
        self.ImposeInletVelocityTemperature()
        debug = True
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ImposeInletVelocityTemperature(self):
        self.outlet_model_part = self.model.GetModelPart(self.outlet_model_part_name)
        for node in self.outlet_model_part.Nodes:
            normal = node.GetValue(KratosMultiphysics.NORMAL)
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            if np.dot(velocity, normal)<0:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.ambient_temperature)
                node.Fix(KratosMultiphysics.TEMPERATURE)
                # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
                # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            debug = True

    def InterpolateVelocityWithMA(self):
        #Import Origin Model Part (Fluid Model Part)
        self.origin_model_part = self.model.CreateModelPart("OriginModelPart")
        self.origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        model_part_io = KratosMultiphysics.ModelPartIO(self.origin_model_part_file_name)
        model_part_io.ReadModelPart(self.origin_model_part) 

        #Set velocity field to origin model part
<<<<<<< Updated upstream
        self.alpha = 1.0
        self.expected_alpha = 1.0
        velocity = self.alpha*np.load("velocity_field.npy")
=======
        velocity = 1.0e-2*np.load("velocity_field.npy")
>>>>>>> Stashed changes
        
        counter = 0
        for node in self.origin_model_part.Nodes:
            node.SetValue(KratosMultiphysics.VELOCITY_X, velocity[counter,0])
            node.SetValue(KratosMultiphysics.VELOCITY_Y, velocity[counter,1])
            node.SetValue(KratosMultiphysics.VELOCITY_Z, velocity[counter,2])
            counter += 1

        #Interpolate velocity to destination model part (Convection-Diffusion)
        # self.destination_model_part = self.model.GetModelPart("ThermalModelPart")
        interpolation = KratosMA.NodalValuesInterpolationProcess2D(self.origin_model_part,self.destination_model_part)
        interpolation.Execute()
        for node in self.destination_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,node.GetValue(KratosMultiphysics.VELOCITY_X))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,node.GetValue(KratosMultiphysics.VELOCITY_Y))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,node.GetValue(KratosMultiphysics.VELOCITY_Z))
        
        debug=True

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
<<<<<<< Updated upstream
        # self.old_alpha = self.alpha
        # self.alpha += 1e-4
        # if self.alpha<self.expected_alpha:
        #     for node in self.destination_model_part.Nodes:
        #         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)/self.old_alpha)
        #         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)/self.old_alpha)
        #         node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)/self.old_alpha)
=======
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
>>>>>>> Stashed changes
        pass

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