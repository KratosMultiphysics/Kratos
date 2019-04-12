import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_lift_process import ComputeLiftProcess

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return ForcesElementsToNodesProcess(Model, settings["Parameters"])

# all the processes python processes should be derived from "python_process"

class ForcesElementsToNodesProcess(ComputeLiftProcess):
    def __init__(self, Model, settings):
        super(ForcesElementsToNodesProcess, self).__init__(Model, settings)

    def ExecuteFinalizeSolutionStep(self):
        self.ComputeForcesOnNodes()
    
    def ComputeForcesOnNodes(self):
        print('COMPUTE LIFT AT NODES')

        KratosMultiphysics.VariableUtils().SetToZero_VectorVar(KratosMultiphysics.FORCE, self.body_model_part.Nodes)

        for cond in self.body_model_part.Conditions:
            n = cond.GetValue(KratosMultiphysics.NORMAL)
            cp = cond.GetValue(KratosMultiphysics.PRESSURE)

            for node in cond.GetNodes():
                added_force = KratosMultiphysics.Vector(3)
                force = KratosMultiphysics.Vector(3)

                added_force = n*(cp/2.0)
                force = node.GetValue(KratosMultiphysics.FORCE) + added_force
                node.SetValue(KratosMultiphysics.FORCE, force)

        total_force = KratosMultiphysics.VariableUtils().SumNonHistoricalNodeVectorVariable(KratosMultiphysics.FORCE, self.body_model_part)
        
        print("Cl = ", total_force[1])
        print("Cd = ", total_force[0])
        print("RZ = ", total_force[2])

        if self.create_output_file:
            with open("cl_points_with_lift.dat", 'w') as cl_file:
                cl_file.write('{0:15.12f}'.format(total_force[1]))


