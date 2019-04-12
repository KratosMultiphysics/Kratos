import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math

def DotProduct(A,B):
    result = 0
    for i in range(len(A)):
        result += A[i]*B[i]
    return result

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
            "reference_area": 1.0,
            "moment_reference_point" : [0.0,0.0,0.0],
            "create_output_file": false
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        self.reference_area =  settings["reference_area"].GetDouble()
        self.moment_reference_point = KratosMultiphysics.Vector(3)
        self.moment_reference_point = settings["moment_reference_point"].GetVector()
        self.create_output_file = settings["create_output_file"].GetBool()
        self.fluid_model_part = self.body_model_part.GetRootModelPart()

    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','COMPUTE LIFT')

        f = KratosMultiphysics.Vector(3)
        m = KratosMultiphysics.Vector(3)

        for cond in self.body_model_part.Conditions:
            n = cond.GetValue(KratosMultiphysics.NORMAL)
            cp = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)

            # Computing forces
            f[0] += n[0]*cp
            f[1] += n[1]*cp
            f[2] += n[2]*cp

            # Computing moment
            mid_point = cond.GetGeometry().Center()
            lever = mid_point-self.moment_reference_point
            m += CrossProduct(lever, n*(-cp))

        Cf = KratosMultiphysics.Vector(3)
        Cf[0] = f[0]/self.reference_area
        Cf[1] = f[1]/self.reference_area
        Cf[2] = f[2]/self.reference_area

        self.Cm = m[2]/self.reference_area

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

        mach_inf = self.fluid_model_part.GetProperties()[1].GetValue(CPFApp.MACH_INFINITY)
        a_inf = self.fluid_model_part.GetProperties()[1].GetValue(KratosMultiphysics.SOUND_VELOCITY)
        u_inf = mach_inf * a_inf
        
        node_velocity_potential_te = te.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
        node_auxiliary_velocity_potential_te = te.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
        if(te.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0):
            potential_jump_phi_minus_psi_te = node_velocity_potential_te - node_auxiliary_velocity_potential_te
        else:
            potential_jump_phi_minus_psi_te = node_auxiliary_velocity_potential_te - node_velocity_potential_te
        self.Cl_te = 2*potential_jump_phi_minus_psi_te/u_inf

    def __ReadWakeDirection(self):
        self.wake_direction = self.fluid_model_part.GetProperties()[1].GetValue(CPFApp.VELOCITY_INFINITY)
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