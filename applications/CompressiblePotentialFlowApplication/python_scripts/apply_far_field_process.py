import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math

def DotProduct(A,B):
    result = 0
    for i in range(len(A)):
        result += A[i]*B[i]
    return result

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "inlet_phi": 1.0,
                "mach_infinity": 0.01,
                "angle_of_attack": 0.0,
                "density_infinity"  : 1.0,
                "speed_of_sound": 340,
                "gamma": 1.4,
                "pressure_infinity": 101325,
                "initialize_flow_field": true
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.far_field_model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.far_field_model_part.GetRootModelPart()

        self.inlet_phi_0 = settings["inlet_phi"].GetDouble()
        self.density_inf = settings["density_infinity"].GetDouble()
        self.mach_inf = settings["mach_infinity"].GetDouble()
        self.gamma = settings["gamma"].GetDouble()
        self.aoa = settings["angle_of_attack"].GetDouble()
        self.pressure_inf = settings["pressure_infinity"].GetDouble()
        self.a_inf = settings["speed_of_sound"].GetDouble()
        self.initialize = settings["initialize_flow_field"].GetBool()

        # Computing free stream velocity
        self.u_inf = self.mach_inf * self.a_inf
        self.velocity_inf = KratosMultiphysics.Vector(3)
        self.velocity_inf[0] = round(self.u_inf*math.cos(self.aoa),8)
        self.velocity_inf[1] = round(self.u_inf*math.sin(self.aoa),8)
        self.velocity_inf[2] = 0.0

        # For the model part
        self.far_field_model_part.ProcessInfo.SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_inf)

        # For the conditions
        self.fluid_model_part.GetProperties()[0].SetValue(CPFApp.DENSITY_INFINITY, self.density_inf)

        # For the elements
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_inf)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.DENSITY_INFINITY, self.density_inf)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.MACH_INFINITY, self.mach_inf)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.GAMMA, self.gamma)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.AOA, self.aoa)
        self.fluid_model_part.GetProperties()[1].SetValue(KratosMultiphysics.SOUND_VELOCITY, self.a_inf)
        self.fluid_model_part.GetProperties()[1].SetValue(CPFApp.PRESSURE_INFINITY, self.pressure_inf)

    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CPFApp.VELOCITY_INFINITY, self.velocity_inf, self.far_field_model_part.Conditions)
        for cond in self.far_field_model_part.Conditions:
            cond.SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_inf)

        # Select the first node in the mesh as reference
        for node in self.far_field_model_part.Nodes:
            x0 = node.X
            y0 = node.Y
            z0 = node.Z
            break

        # Find smallest distance_to_reference
        pos = 1e30
        for node in self.far_field_model_part.Nodes:
            # Computing distance to reference
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            distance_to_reference = dx*self.velocity_inf[0] + dy*self.velocity_inf[1] + dz*self.velocity_inf[2]

            if(distance_to_reference < pos):
                pos = distance_to_reference
                self.reference_inlet_node = node

        # Fix nodes in the inlet
        for cond in self.far_field_model_part.Conditions:
            normal = cond.GetValue(KratosMultiphysics.NORMAL)
            projection = DotProduct(normal,self.velocity_inf)

            if( projection < 0):
                for node in cond.GetNodes():
                    # Computing distance to reference
                    dx = node.X - self.reference_inlet_node.X
                    dy = node.Y - self.reference_inlet_node.Y
                    dz = node.Z - self.reference_inlet_node.Z

                    inlet_phi = dx*self.velocity_inf[0] + dy*self.velocity_inf[1] + dz*self.velocity_inf[2]
                    node.Fix(CPFApp.VELOCITY_POTENTIAL)
                    node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,inlet_phi + self.inlet_phi_0)
                    if self.far_field_model_part.HasNodalSolutionStepVariable(CPFApp.ADJOINT_VELOCITY_POTENTIAL):
                        node.Fix(CPFApp.ADJOINT_VELOCITY_POTENTIAL)
                        node.SetSolutionStepValue(CPFApp.ADJOINT_VELOCITY_POTENTIAL,0,inlet_phi)

        if(self.initialize):
            for node in self.fluid_model_part.Nodes:
                # Computing distance to reference
                dx = node.X - self.reference_inlet_node.X
                dy = node.Y - self.reference_inlet_node.Y
                dz = node.Z - self.reference_inlet_node.Z

                initial_phi = dx*self.velocity_inf[0] + dy*self.velocity_inf[1] + dz*self.velocity_inf[2]
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,initial_phi + self.inlet_phi_0)
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL,0,initial_phi + self.inlet_phi_0)
                #TODO: How to initialize the adjoint potential field?
                '''
                if self.far_field_model_part.HasNodalSolutionStepVariable(CPFApp.ADJOINT_VELOCITY_POTENTIAL):
                        node.SetSolutionStepValue(CPFApp.ADJOINT_VELOCITY_POTENTIAL,0,initial_phi)
                        node.SetSolutionStepValue(CPFApp.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0,initial_phi)
                '''

    def ExecuteInitializeSolutionStep(self):
        self.Execute()
