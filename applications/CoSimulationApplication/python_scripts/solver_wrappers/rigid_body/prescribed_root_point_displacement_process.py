
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.prescribed_value_process import _PrescribedValueProcess
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC


def Factory(solver, parameters):
    return PrescribedRootPointDisplacementProcess(solver, parameters)


class PrescribedRootPointDisplacementProcess(_PrescribedValueProcess):
    
    def _ApplyValue(self, index, prescribed_value):
        if index < 3:
            variable = KMC.PRESCRIBED_DISPLACEMENT
        else:
            index -= 3
            variable = KMC.PRESCRIBED_ROTATION
        values = self.solver.root_point_model_part.Nodes[2].GetSolutionStepValue(variable, 0)
        values[index] += prescribed_value
        self.solver.root_point_model_part.Nodes[2].SetSolutionStepValue(variable, 0, values)
        #self.solver.prescribed_root_point_displ[index] += prescribed_value