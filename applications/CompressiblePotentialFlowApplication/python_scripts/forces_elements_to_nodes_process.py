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
        super(ForcesElementsToNodesProcess, self).ExecuteFinalizeSolutionStep()
        self.ComputeForcesOnNodes()
    
    def ComputeForcesOnNodes(self):
        initialize_force = KratosMultiphysics.Vector(3)
        for node in self.body_model_part.Nodes:
            node.SetValue(KratosMultiphysics.FORCE, initialize_force)

        for cond in self.body_model_part.Conditions:
            n = cond.GetValue(KratosMultiphysics.NORMAL)
            cp = cond.GetValue(KratosMultiphysics.PRESSURE)

            for node in cond.GetNodes():
                tmp = node.GetValue(KratosMultiphysics.FORCE)
                added_force = KratosMultiphysics.Vector(3)
                added_force[0] = n[0]*cp/2.0
                added_force[1] = n[1]*cp/2.0
                added_force[2] = n[2]*cp/2.0
                force = KratosMultiphysics.Vector(3)
                force = tmp + added_force
                node.SetValue(KratosMultiphysics.FORCE, force)

        total_force = KratosMultiphysics.Vector(3)
        for node in self.body_model_part.Nodes:
            force = node.GetValue(KratosMultiphysics.FORCE)
            total_force += force

        if self.create_output_file:
            with open("cl_points_with_lift.dat", 'w') as cl_file:
                cl_file.write('{0:15.12f}'.format(total_force[1]))

        print("Forces Transfered To Nodes")

