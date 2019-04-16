import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math

def DotProduct(A,B):
    return sum(i[0]*i[1] for i in zip(A, B))

def CrossProduct(A, B):
    C = KratosMultiphysics.Vector(3)
    C[0] = A[1]*B[2]-A[2]*B[1]
    C[1] = A[2]*B[0]-A[0]*B[2]
    C[2] = A[0]*B[1]-A[1]*B[0]
    return C

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "please specify the model part that contains the surface nodes",
            "velocity_infinity": [1.0,0.0,0.0],
            "reference_area": 1.0,
            "moment_reference_point" : [0.0,0.0,0.0],
            "create_output_file": false
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.body_model_part.GetRootModelPart()
        self.reference_area =  settings["reference_area"].GetDouble()
        self.create_output_file = settings["create_output_file"].GetBool()
        self.moment_reference_point = KratosMultiphysics.Vector(3)
        self.moment_reference_point = settings["moment_reference_point"].GetVector()
        

    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','COMPUTE LIFT')

        total_force = KratosMultiphysics.Vector(3)
        total_moment = KratosMultiphysics.Vector(3)

        for cond in self.body_model_part.Conditions:
            n = cond.GetValue(KratosMultiphysics.NORMAL)
            cp = cond.GetValue(KratosMultiphysics.PRESSURE)

            # Computing forces
            total_force += n.__mul__(cp)

            # Computing moment
            mid_point = cond.GetGeometry().Center()
            lever = mid_point-self.moment_reference_point
            total_moment += CrossProduct(lever, n*(-cp))

        Cf = total_force.__div__(self.reference_area)
        self.Cm = total_moment[2]/self.reference_area

        self.__ReadWakeDirection()

        self.Cl = DotProduct(Cf,self.wake_normal)
        self.Cd = DotProduct(Cf,self.wake_direction)

        self.__ComputeLiftJump()

        KratosMultiphysics.Logger.PrintInfo(' Cl = ', self.Cl)
        KratosMultiphysics.Logger.PrintInfo(' Cd = ', self.Cd)
        KratosMultiphysics.Logger.PrintInfo(' RZ = ', Cf[2])
        KratosMultiphysics.Logger.PrintInfo(' Cm = ', self.Cm)
        KratosMultiphysics.Logger.PrintInfo(' Cl = ' , self.Cl_te, ' = 2 * DPhi / U_inf ')

        if self.create_output_file:
            with open("cl.dat", 'w') as cl_file:
                cl_file.write('{0:15.12f}'.format(self.Cl))
            with open("moment.dat", 'w') as mom_file:
                mom_file.write('{0:15.12f}'.format(self.Cm))
                with open("cl_jump.dat", 'w') as cl_file:
                 cl_file.write('{0:15.12f}'.format(self.Cl_te))

    def __ComputeLiftJump(self):
        # Find the Trailing Edge node
        for node in self.body_model_part.Nodes:
            if node.GetValue(CPFApp.TRAILING_EDGE):
                 te=node
                 break

        velocity_infinity = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.VELOCITY_INFINITY)
        u_inf = math.sqrt(DotProduct(velocity_infinity,velocity_infinity))
        
        node_velocity_potential_te = te.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
        node_auxiliary_velocity_potential_te = te.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
        if(te.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0):
            potential_jump_phi_minus_psi_te = node_velocity_potential_te - node_auxiliary_velocity_potential_te
        else:
            potential_jump_phi_minus_psi_te = node_auxiliary_velocity_potential_te - node_velocity_potential_te
        self.Cl_te = 2*potential_jump_phi_minus_psi_te/u_inf

    def __ReadWakeDirection(self):
        self.wake_direction = self.fluid_model_part.ProcessInfo.GetValue(CPFApp.VELOCITY_INFINITY)
        if(self.wake_direction.Size() != 3):
            raise Exception('The wake direction should be a vector with 3 components!')

        dnorm = math.sqrt(
            self.wake_direction[0]**2 + self.wake_direction[1]**2 + self.wake_direction[2]**2)
        self.wake_direction[0] /= dnorm
        self.wake_direction[1] /= dnorm
        self.wake_direction[2] /= dnorm

        self.wake_normal = KratosMultiphysics.Vector(3)
        self.wake_normal[0] = -self.wake_direction[1]
        self.wake_normal[1] = self.wake_direction[0]
        self.wake_normal[2] = 0.0