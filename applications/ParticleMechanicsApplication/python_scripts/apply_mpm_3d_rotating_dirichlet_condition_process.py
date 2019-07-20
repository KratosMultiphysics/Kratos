import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
from KratosMultiphysics.ParticleMechanicsApplication.apply_mpm_particle_dirichlet_condition_process import ApplyMPMParticleDirichletConditionProcess

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
                "particles_per_condition"   : 0,
                "imposition_type"           : "penalty",
                "penalty_factor"            : 0,
                "constrained"               : "fixed",
                "option"                    : "",
                "rotation_center"           : [0.0, 0.0, 0.0],
                "rotation_velocity"         : [0.0, 0.0, 0.0],
                "compute_rotation_center"   : false,
                "rotation_option"           : ""
            }  """ )

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
        super(ApplyMPM3DRotatingDirichletConditionProcess, self).__init__(Model, settings)


    def ExecuteBeforeSolutionLoop(self):
        # Get updated model_part
        if (self.model_part_name.startswith('Background_Grid.')):
            self.model_part_name = self.model_part_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.model_part_name
        self.model_part = self.model[mpm_material_model_part_name]

        # If compute rotation center automatically is needed
        if self.compute_rotation_center:
            self._ComputeCenterRotation()

        # Compute initial radius
        self._ComputeRadius()

        # Initial Quaternion, its increment, and rotation_matrix
        self._InitializeRotationVariables()

    def ExecuteInitializeSolutionStep(self):

        # If compute rotation center automatically is needed and if the center is moved
        if (self.compute_rotation_center and self.rotation_option == "moving_center"):
            self._ComputeCenterRotation()

        # Compute delta quaternion
        self._ComputeDeltaQuaternion(self.rotation_velocity)

        # Compute quaternion multiplication
        self._quaternion = self._MultiplyQuaternion(self._delta_quaternion, self._quaternion)

        # Normalize quaternion
        self._NormalizeQuaternion(self._quaternion)

        # Compute rotation matrix
        new_rotation_matrix = self._ComputeRotationMatrixFromQuaternion(self._quaternion)

        for mpc in self.model_part.Conditions:
            # Update impose_displacement
            imposed_disp = mpc.GetValue(KratosParticle.MPC_IMPOSED_DISPLACEMENT)
            initial_radius = mpc.GetValue(KratosParticle.MPC_RADIUS)
            imposed_disp += new_rotation_matrix * initial_radius - self._rotation_matrix * initial_radius
            mpc.SetValue(KratosParticle.MPC_IMPOSED_DISPLACEMENT, imposed_disp)

            # Update normal vector
            normal = mpc.GetValue(KratosParticle.MPC_NORMAL)
            modified_normal = new_rotation_matrix * (self._TransposeMatrix(self._rotation_matrix) * normal)
            mpc.SetValue(KratosParticle.MPC_NORMAL, modified_normal)

        # Copy rotation matrix
        self._rotation_matrix = new_rotation_matrix

    ### Protected functions
    def _ComputeCenterRotation(self):
        auto_rc = KratosMultiphysics.Vector(3)
        for mpc in self.model_part.Conditions:
            auto_rc += mpc.GetValue(KratosParticle.MPC_COORD)

        auto_rc = auto_rc / self.model_part.NumberOfConditions()
        self.rotation_center = auto_rc

    def _ComputeRadius(self):
        for mpc in self.model_part.Conditions:
            xg_c = mpc.GetValue(KratosParticle.MPC_COORD)
            radius = xg_c - self.rotation_center
            mpc.SetValue(KratosParticle.MPC_RADIUS, radius)

    def _InitializeRotationVariables(self):
        self._quaternion = KratosMultiphysics.Vector(4)
        self._quaternion[0] = 1.0
        self._quaternion[1] = 0.0
        self._quaternion[2] = 0.0
        self._quaternion[3] = 0.0

        self._delta_quaternion = self._quaternion

        self._rotation_matrix = KratosMultiphysics.Matrix(3,3)
        self._rotation_matrix.fill(0.0)
        self._rotation_matrix[0,0] = 1.0
        self._rotation_matrix[1,1] = 1.0
        self._rotation_matrix[2,2] = 1.0

    def _ComputeDeltaQuaternion(self, rot_velocity):
        if((not isinstance(rot_velocity, KratosMultiphysics.Vector)) and rot_velocity.Size() != 3):
            raise Exception("expected input shall be a KratosMultiphysics.Vector object with size 3")

        omega = math.sqrt(rot_velocity[0]*rot_velocity[0] + rot_velocity[1]*rot_velocity[1] + rot_velocity[2]*rot_velocity[2])
        if omega < 1.e-12:
            omega = 1.e-12

        dt = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        self._delta_quaternion[0] = math.cos(omega * dt / 2.0)
        self._delta_quaternion[1] = math.sin(omega * dt / 2.0) * rot_velocity[0] / omega
        self._delta_quaternion[2] = math.sin(omega * dt / 2.0) * rot_velocity[1] / omega
        self._delta_quaternion[3] = math.sin(omega * dt / 2.0) * rot_velocity[2] / omega

    def _MultiplyQuaternion(self, p, q):
        if((not isinstance(p, KratosMultiphysics.Vector)) and p.Size() != 4):
            raise Exception("expected input: p shall be a KratosMultiphysics.Vector object with size 4")
        if((not isinstance(q, KratosMultiphysics.Vector)) and q.Size() != 4):
            raise Exception("expected input: q shall be a KratosMultiphysics.Vector object with size 4")

        result = KratosMultiphysics.Vector(4)
        result[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]
        result[1] = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2]
        result[2] = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1]
        result[3] = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0]

        return result

    def _NormalizeQuaternion(self, q):
        if((not isinstance(q, KratosMultiphysics.Vector)) and q.Size() != 4):
            raise Exception("expected input shall be a KratosMultiphysics.Vector object with size 4")

        len_q = math.sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])
        if (len_q < 1.e-16):
            len_q = 1.e-16

        q = q/len_q

    def _ComputeRotationMatrixFromQuaternion(self, q):
        if((not isinstance(q, KratosMultiphysics.Vector)) and q.Size() != 4):
            raise Exception("expected input shall be a KratosMultiphysics.Vector object with size 4")

        rot_mat = KratosMultiphysics.Matrix(3,3)
        rot_mat[0,0] = 2.0 * (q[0]*q[0] + q[1]*q[1]) - 1.0
        rot_mat[0,1] = 2.0 * (q[1]*q[2] - q[0]*q[3])
        rot_mat[0,2] = 2.0 * (q[1]*q[3] + q[0]*q[2])

        rot_mat[1,0] = 2.0 * (q[1]*q[2] + q[0]*q[3])
        rot_mat[1,1] = 2.0 * (q[0]*q[0] + q[2]*q[2]) - 1.0
        rot_mat[1,2] = 2.0 * (q[2]*q[3] - q[0]*q[1])

        rot_mat[2,0] = 2.0 * (q[1]*q[3] - q[0]*q[2])
        rot_mat[2,1] = 2.0 * (q[2]*q[3] + q[0]*q[1])
        rot_mat[2,2] = 2.0 * (q[0]*q[0] + q[3]*q[3]) - 1.0

        return rot_mat

    def _TransposeMatrix(self, A):
        if(not isinstance(A, KratosMultiphysics.Matrix)):
            raise Exception("expected input shall be a KratosMultiphysics.Matrix object with equal number of row-column (square)")

        result = KratosMultiphysics.Matrix(A.Size2(),A.Size1())
        for i in range (A.Size1()):
            for j in range (A.Size2()):
                result[j,i] = A[i,j]

        return result