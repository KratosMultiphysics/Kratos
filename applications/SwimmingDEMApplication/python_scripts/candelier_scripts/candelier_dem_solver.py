from KratosMultiphysics import *
import swimming_DEM_solver
import candelier_scripts.candelier_parameters as candelier_pp

BaseSolver = swimming_DEM_solver.SwimmingDEMSolver

def Cross(a, b):
    c0 = a[1]*b[2] - a[2]*b[1]
    c1 = a[2]*b[0] - a[0]*b[2]
    c2 = a[0]*b[1] - a[1]*b[0]
    return Vector([c0, c1, c2])

class CandelierDEMSolver(BaseSolver):
    def __init__(self, model, project_parameters, fluid_solver, dem_solver, pp):
        super(CandelierDEMSolver, self).__init__(model, project_parameters, fluid_solver, dem_solver, pp)
        self.frame_angular_vel = Vector([0, 0, self.project_parameters["angular_velocity_of_frame_Z"].GetDouble()])
        self.omega = self.project_parameters["angular_velocity_of_frame_Z"].GetDouble()

    def SolveDEMSolutionStep(self):
        super(CandelierDEMSolver, self).SolveDEMSolutionStep()

    def ApplyForwardCoupling(self, alpha = 'None'):
        super(CandelierDEMSolver, self).ApplyForwardCoupling(alpha)

        for node in self.dem_solver.spheres_model_part.Nodes:
            omega = candelier_pp.omega
            r = Vector([node.X, node.Y, node.Z])
            vx = - omega * r[1]
            vy =   omega * r[0]
            v = Vector([vx, vy, 0])
            ax = - r[0] * omega ** 2
            ay = - r[1] * omega ** 2
            a = Vector([ax, ay, 0])

            if self.is_rotating_frame:
                v = self.GetVelocityRelativeToMovingFrame(r_rel = r, v_glob = v)
                a = self.GetAccelerationRelativeToMovingFrame(r_rel = r, v_rel = v, a_glob = a)

            node.SetSolutionStepValue(FLUID_VEL_PROJECTED, v)
            node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED, a)

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self, time=None):
        super(CandelierDEMSolver, self).ApplyForwardCouplingOfVelocityToSlipVelocityOnly()
        for node in self.dem_solver.spheres_model_part.Nodes:
            r = Vector([node.X, node.Y, node.Z])
            new_vx = - candelier_pp.omega * r[1]
            new_vy =   candelier_pp.omega * r[0]
            new_v = Vector([new_vx, new_vy, 0.])

            if self.is_rotating_frame:
                new_v = self.GetVelocityRelativeToMovingFrame(r_rel = r, v_glob = new_v)

            # the current FLUID_VEL_PROJECTED is still needed and so we use
            # SLIP_VELOCITY to store it.
            node.SetSolutionStepValue(SLIP_VELOCITY, new_v)

    def GetVelocityRelativeToMovingFrame(self, r_rel, v_glob):
        cross_omega_r = Cross(self.frame_angular_vel, r_rel)
        return v_glob - cross_omega_r

    def GetAccelerationRelativeToMovingFrame(self, r_rel, v_rel, a_glob):
        cross_omega_r = Cross(self.frame_angular_vel, r_rel)
        cross_omega_v = Cross(self.frame_angular_vel, v_rel)
        cross_omega_omega_r = Cross(self.frame_angular_vel, cross_omega_r)
        return a_glob - 2 * cross_omega_v - cross_omega_omega_r

    def ImportModelPart(self): # TODO: implement this
        pass