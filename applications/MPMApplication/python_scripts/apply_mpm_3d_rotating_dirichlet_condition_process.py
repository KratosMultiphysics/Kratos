import KratosMultiphysics
from KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.apply_mpm_particle_dirichlet_condition_process import ApplyMPMParticleDirichletConditionProcess

import math

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPM3DRotatingDirichletConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPM3DRotatingDirichletConditionProcess(ApplyMPMParticleDirichletConditionProcess):
    def __init__(self, Model, settings ):

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "material_points_per_condition"   : 0,
                "imposition_type"           : "penalty",
                "penalty_factor"            : 0,
                "constrained"               : "fixed",
                "option"                    : "",
                "rotation_center"           : [0.0, 0.0, 0.0],
                "rotation_velocity"         : [0.0, 0.0, 0.0],
                "compute_rotation_center"   : false,
                "rotation_option"           : ""
            }  """ )

        context_string = type(self).__name__
        old_name = 'particles_per_condition'
        new_name = 'material_points_per_condition'
        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.model = Model
        self.model_part_name = settings["model_part_name"].GetString()

        # Private process variables
        self.rotation_center = settings["rotation_center"].GetVector()
        self.rotation_velocity = settings["rotation_velocity"].GetVector()
        self.compute_rotation_center = settings["compute_rotation_center"].GetBool()
        self.rotation_option = settings["rotation_option"].GetString()

        settings.RemoveValue("rotation_center")
        settings.RemoveValue("rotation_velocity")
        settings.RemoveValue("compute_rotation_center")
        settings.RemoveValue("rotation_option")

        # Initiate base class - Dirichlet condition
        super().__init__(Model, settings)


    def ExecuteBeforeSolutionLoop(self):
        # Get updated model_part
        if (self.model_part_name.startswith('Background_Grid.')):
            self.model_part_name = self.model_part_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.model_part_name
        self.model_part = self.model[mpm_material_model_part_name]

        # If compute rotation center automatically is needed
        if self.compute_rotation_center:
            self._ComputeCenterRotation()

        # Initial Quaternion, its increment, and rotation_matrix
        self._InitializeRotationVariables()

    def ExecuteInitializeSolutionStep(self):

        # If compute rotation center automatically is needed and if the center is moved
        if self.compute_rotation_center and self.rotation_option == "moving_center":
            self._ComputeCenterRotation()

        # Compute delta quaternion
        self._ComputeDeltaQuaternion(self.rotation_velocity)
        
        self._quaternion = self._quaternion.MultiplyQuaternion(self._delta_quaternion, self._quaternion)
        
        # Normalize quaternion
        self._quaternion.Normalize()
        
        # Compute rotation matrix
        new_rotation_matrix = KratosMultiphysics.Matrix(3,3)
        self._quaternion.ToRotationMatrix(new_rotation_matrix)

        for mpc in self.model_part.Conditions:
            # Compute current radius
            mpc_coord = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD,self.model_part.ProcessInfo)[0]
            radius = mpc_coord - self.rotation_center

            # Update impose_displacement
            imposed_disp = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_IMPOSED_DISPLACEMENT,self.model_part.ProcessInfo)[0]
            
            imposed_disp += new_rotation_matrix * (self._rotation_matrix.transpose() * radius) - radius
            mpc.SetValuesOnIntegrationPoints(KratosMPM.MPC_IMPOSED_DISPLACEMENT,[imposed_disp],self.model_part.ProcessInfo)

            # Update normal vector
            normal = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_NORMAL,self.model_part.ProcessInfo)[0]
            modified_normal = new_rotation_matrix * (self._rotation_matrix.transpose() * normal)
            mpc.SetValuesOnIntegrationPoints(KratosMPM.MPC_NORMAL,[modified_normal],self.model_part.ProcessInfo)


        # Copy rotation matrix
        self._rotation_matrix = new_rotation_matrix

    ### Protected functions
    def _ComputeCenterRotation(self):
        auto_rc = KratosMultiphysics.Vector(3)
        for mpc in self.model_part.Conditions:
            auto_rc += mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD,self.model_part.ProcessInfo)[0]

        auto_rc = auto_rc / self.model_part.NumberOfConditions()
        self.rotation_center = auto_rc

    def _InitializeRotationVariables(self):
        self._quaternion = KratosMultiphysics.Quaternion()
        self._quaternion = self._quaternion.Identity()
        
        self._delta_quaternion = KratosMultiphysics.Quaternion()
        self._delta_quaternion = self._delta_quaternion.Identity()

        self._rotation_matrix = KratosMultiphysics.Matrix(3,3)
        self._rotation_matrix.fill_identity()

    def _ComputeDeltaQuaternion(self, rot_velocity):
        if (not isinstance(rot_velocity, KratosMultiphysics.Vector)) and rot_velocity.Size() != 3:
            raise Exception("expected input shall be a KratosMultiphysics.Vector object with size 3")

        omega = math.sqrt(rot_velocity[0]*rot_velocity[0] + rot_velocity[1]*rot_velocity[1] + rot_velocity[2]*rot_velocity[2])
        if omega < 1.e-12:
            omega = 1.e-12

        dt = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        
        self._delta_quaternion = KratosMultiphysics.Quaternion()
        self._delta_quaternion.W = math.cos(omega * dt / 2.0)
        self._delta_quaternion.X = math.sin(omega * dt / 2.0) * rot_velocity[0] / omega
        self._delta_quaternion.Y = math.sin(omega * dt / 2.0) * rot_velocity[1] / omega
        self._delta_quaternion.Z = math.sin(omega * dt / 2.0) * rot_velocity[2] / omega
