import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D as DefineWakeProcess2D

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return ComputeLiftJumpProcess2D(Model, settings["Parameters"])

# TODO Implement this process in C++ and make it open mp parallel to save time selecting the wake elements
class ComputeLiftJumpProcess2D(DefineWakeProcess2D):
    def __init__(self, Model, settings):
        # Call the base Kratos process constructor
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "velocity_infinity": [1.0,0.0,0],
            "create_output_file": false
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        body_model_part_name = settings["model_part_name"].GetString()
        if body_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess2D\n"
            err_msg += "Please specify the model part that contains the body surface nodes"
            raise Exception(err_msg)
        self.body_model_part = Model[body_model_part_name]

        self.velocity_infinity = [0,0,0]
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.create_output_file = settings["create_output_file"].GetBool()

    def ExecuteInitialize(self):
        # Save the trailing edge for further computations
        self.SaveTrailingEdgeNode()

    def ExecuteFinalizeSolutionStep(self):
        node_velocity_potential_te = self.te.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
        node_auxiliary_velocity_potential_te = self.te.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
        if(self.te.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0):
            potential_jump_phi_minus_psi_te = node_velocity_potential_te - node_auxiliary_velocity_potential_te
        else:
            potential_jump_phi_minus_psi_te = node_auxiliary_velocity_potential_te - node_velocity_potential_te
        Cl_te = 2*potential_jump_phi_minus_psi_te/self.velocity_infinity[0]
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftJumpProcess2D','Potential Jump: Phi - Psi (at trailing edge node) = ', potential_jump_phi_minus_psi_te, '=> CL = ',Cl_te)

        if self.create_output_file:
             with open("cl_jump.dat", 'w') as cl_file:
                 cl_file.write('{0:15.12f}'.format(Cl_te))
