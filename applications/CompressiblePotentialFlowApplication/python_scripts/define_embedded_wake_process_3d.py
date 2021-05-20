import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import time
import math


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineEmbeddedWakeProcess(Model, settings["Parameters"])


class DefineEmbeddedWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        # Call the base Kratos process constructor
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "epsilon": 1e-9
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.main_model_part = Model[settings["model_part_name"].GetString()].GetRootModelPart()
        self.wake_model_part = Model.CreateModelPart("wake")
        # self.plane_model_part = Model.CreateModelPart("plane")
        # self.output_model_part = Model.CreateModelPart("output")

        self.epsilon = settings["epsilon"].GetDouble()

    def ExecuteInitialize(self):
        ini_time = time.time()

        self._DefineWakeModelPart()

        # self._MoveAndRotateWake()
        # Executing define wake process
        # KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess3D(self.main_model_part, self.wake_model_part).Execute()
        CPFApp.DefineEmbeddedWakeProcess3D(self.main_model_part, self.wake_model_part).Execute()

        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Wake computation time: ',time.time()-ini_time)

    def _DefineWakeModelPart(self):
        ''' This function generates the modelpart of the wake. TODO: make end of the domain user-definable.
        '''
        self.wake_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.wake_model_part.CreateNewNode(2, 0.0, 2.0, 0.0)

        self.wake_model_part.CreateNewNode(3, 200.0*math.cos(math.radians(-5.0)), 0.5, 200.0*math.sin(math.radians(-5.0)))
        self.wake_model_part.CreateNewNode(4, 200.0*math.cos(math.radians(-5.0)), 0.0, 200.0*math.sin(math.radians(-5.0)))

        # self.wake_model_part.CreateNewNode(1, 0.4980874, -0.1, -0.04357787)
        # self.wake_model_part.CreateNewNode(2, 0.4980874, 0.6, -0.04357787)
        # self.wake_model_part.CreateNewNode(3, 200.0, 0.6, -0.04357787)
        # self.wake_model_part.CreateNewNode(4, 200.0, -0.1, -0.04357787)

        self.wake_model_part.CreateNewElement("Element3D3N", 1, [1,3,2], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewElement("Element3D3N", 2, [1,4,3], KratosMultiphysics.Properties(0))

    def _MoveAndRotateWake(self):
        ''' This function moves and rotates the wake with the same parameters as the geometry.
        '''
        self.moving_parameters = KratosMultiphysics.Parameters()
        self.moving_parameters.AddEmptyValue("origin")
        self.moving_parameters["origin"].SetVector(self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        self.moving_parameters.AddEmptyValue("rotation_angle")
        angle=math.radians(-self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE))
        self.moving_parameters["rotation_angle"].SetDouble(angle)
        CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()

    # def ExecuteFinalizeSolutionStep(self):
    #     CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(self.wake_sub_model_part, 1e-1, 0)