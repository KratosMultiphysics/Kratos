import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return ComputeForcesOnNodesProcess(Model, settings["Parameters"])

# all the processes python processes should be derived from "python_process"

class ComputeForcesOnNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "please specify the model part that contains the surface nodes",
            "create_output_file": false
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        self.create_output_file = settings["create_output_file"].GetBool()

    def ExecuteFinalizeSolutionStep(self):
        self.Execute()

    def Execute(self):
        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess', 'COMPUTE REACTIONS ON NODES')

        # density_infinity = self.body_model_part.GetProperties()[0].Has(DENSITY)
        # KratosMultiphysics.Logger.PrintInfo('DENSITY INFINITY', density_infinity)

        KratosMultiphysics.VariableUtils().SetToZero_VectorVar(KratosMultiphysics.REACTION, self.body_model_part.Nodes)

        velocity_infinity = self.body_model_part.ProcessInfo.GetValue(CPFApp.VELOCITY_INFINITY)
        self.body_model_part.ProcessInfo.SetValue(KratosMultiphysics.DENSITY, 1.225)
        density_infinity = self.body_model_part.ProcessInfo.GetValue(KratosMultiphysics.DENSITY)
        u_inf = velocity_infinity.norm_2()
        dynamic_pressure = 0.5*density_infinity*u_inf**2

        for cond in self.body_model_part.Conditions:
            n = cond.GetValue(KratosMultiphysics.NORMAL)
            cp = cond.GetValue(KratosMultiphysics.PRESSURE)

            for node in cond.GetNodes():
                added_force = n*(cp/2.0)*dynamic_pressure
                force = node.GetValue(KratosMultiphysics.REACTION) + added_force
                node.SetValue(KratosMultiphysics.REACTION, force)

        total_force = KratosMultiphysics.VariableUtils().SumNonHistoricalNodeVectorVariable(KratosMultiphysics.REACTION, self.body_model_part)

        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess','Lift Force = ', total_force[1])
        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess','Drag Force = ', total_force[0])
        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess','Side Force = ', total_force[2])

        if self.create_output_file:
            with open("cl_points_with_lift.dat", 'w') as cl_file:
                cl_file.write('{0:15.12f}'.format(total_force[1]/dynamic_pressure))

