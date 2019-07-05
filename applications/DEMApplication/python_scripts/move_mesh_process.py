from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as KratosDem


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MoveMeshProcess(model, settings["Parameters"])

class MoveMeshProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"                : "",
            "time_step"                      : 0.1,
            "previous_displ"                 : [0.0, 0.0, 0.0],
            "linear_velocity"                : [0.0, 0.0, 0.0],
            "velocity_start_time"            : 0.0,
            "velocity_stop_time "            : 1.0,
            "linear_period "                 : 1.0,
            "angular_velocity"               : [0.0, 0.0, 0.0],
            "rotation_center"                : [0.0, 0.0, 0.0],
            "angular_velocity_start_time"    : 0.0,
            "angular_velocity_stop_time "    : 1.0,
            "angular_velocity_period "       : 0.0,
            "rigid_body_motion"              : true,
            "fixed_mesh "                    : false
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[self.params["model_part_name"].GetString()]

        self.time_step = self.params["time_step"].GetDouble()

        self.previous_displ = [0.0,0.0,0.0]
        self.linear_velocity = self.params["linear_velocity"].GetVector()
        self.velocity_start_time = self.params["velocity_start_time"].GetDouble()
        self.angular_velocity = self.params["angular_velocity"].GetVector()
        self.angular_velocity_start_time = self.params["angular_velocity_start_time"].GetDouble()
        self.rigid_body_motion = self.params["rigid_body_motion"].GetBool()
        self.rotation_center = self.params["rotation_center"].GetVector()
        self.fixed_mesh = False #self.params["fixed_mesh"].GetBool()
        self.angular_velocity_stop_time = 0.0#self.params["angular_velocity_stop_time"].GetDouble()
        self.velocity_stop_time = 1000.0 #self.params["velocity_stop_time"].GetDouble()
        self.angular_velocity_period = 0.0 # self.params["angular_velocity_period"].GetDouble()
        self.linear_period = 1.0#self.params["linear_period"].GetDouble()
        self.angular_period = 0.0#self.params["linear_period"].GetDouble()
        self.mod_angular_velocity = 50#sqrt(self.angular_velocity[0]*self.angular_velocity[0] + self.angular_velocity[1]*self.angular_velocity[1] + self.angular_velocity[2]*self.angular_velocity[2])
        self.initial_center = [0,0,0]
    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        center_position = [0,0,0]
        linear_velocity_changed = [0,0,0]
        angular_velocity_changed = [0,0,0]
        new_axes1 = [0,0,0]
        new_axes2 = [0,0,0]
        new_axes3 = [0,0,0]
        
        print(self.model_part.GetNode(4))
        KratosDem.TranslateGridOfNodes(
            time,
            self.velocity_start_time,
            self.velocity_stop_time,
            center_position,
            self.rotation_center,
            self.previous_displ,
            linear_velocity_changed,
            self.linear_period,
            self.time_step,
            self.linear_velocity)

        KratosDem.RotateGridOfNodes(
            time,
            self.angular_velocity_start_time,
            self.angular_velocity_stop_time,
            angular_velocity_changed,
            self.angular_period,
            self.mod_angular_velocity,
            self.angular_velocity,
            new_axes1,
            new_axes2,
            new_axes3)

        KratosDem.UpdateKinematicVariablesOfAGridOfNodes(
             self.mod_angular_velocity,
             self.linear_velocity,
             self.initial_center,
             new_axes1,
             new_axes2,
             new_axes3,
             angular_velocity_changed,
             linear_velocity_changed,
             center_position,
             self.fixed_mesh,
             self.time_step,
             self.model_part.Nodes);


        print(self.model_part.GetNode(4))