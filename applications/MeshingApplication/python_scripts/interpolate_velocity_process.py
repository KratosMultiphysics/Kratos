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
            "destination_model_part_name" : "ModelPartName"
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model = model
        self.settings = settings

    def ExecuteInitialize(self):
        self.origin_model_part_file_name = self.settings["origin_model_part_file_name"].GetString()
        destination_model_part_name = self.settings["destination_model_part_name"].GetString()
        self.destination_model_part = self.model.GetModelPart(destination_model_part_name)
        
        self.InterpolateVelocityWithMA()

    def InterpolateVelocityWithMA(self):
        #Import Origin Model Part 
        self.origin_model_part = self.model.CreateModelPart("OriginModelPart")
        for variable in self.destination_model_part.GetHistoricalVariablesNames():
            print(variable)
            if KratosMultiphysics.KratosGlobals.HasVariable(variable):
                self.origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))

        model_part_io = KratosMultiphysics.ModelPartIO(self.origin_model_part_file_name)
        model_part_io.ReadModelPart(self.origin_model_part) 

        #Set velocity field to origin model part
        self.alpha = 1e-3
        self.expected_alpha = 1.0
        
        for node in self.origin_model_part.Nodes:
            node.SetValue(KratosMultiphysics.VELOCITY_X, self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
            node.SetValue(KratosMultiphysics.VELOCITY_Y, self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
            node.SetValue(KratosMultiphysics.VELOCITY_Z, self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z))
            node.SetValue(KratosMultiphysics.PRESSURE, self.alpha*node.GetSolutionStepValue(KratosMultiphysics.PRESSURE))

        #Interpolate velocity to destination model part 
        interpolation = KratosMA.NodalValuesInterpolationProcess2D(self.origin_model_part,self.destination_model_part)
        interpolation.Execute()
        for node in self.destination_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,node.GetValue(KratosMultiphysics.VELOCITY_X))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,node.GetValue(KratosMultiphysics.VELOCITY_Y))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,node.GetValue(KratosMultiphysics.VELOCITY_Z))
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,node.GetValue(KratosMultiphysics.PRESSURE))

    def ExecuteInitializeSolutionStep(self):
        self.old_alpha = self.alpha
        self.alpha += 1e-3
        if self.alpha<self.expected_alpha:
            for node in self.destination_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)/self.old_alpha)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)/self.old_alpha)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,self.alpha*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)/self.old_alpha)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,self.alpha*node.GetValue(KratosMultiphysics.PRESSURE)/self.old_alpha)
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)
                node.Fix(KratosMultiphysics.PRESSURE)