
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.prescribed_value_process import _PrescribedValueProcess
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC


def Factory(solver, parameters):
    return PrescribedForceProcess(solver, parameters)


class PrescribedForceProcess(_PrescribedValueProcess):
    
    def _ApplyValue(self, index, prescribed_value):
        if index < 3:
            variable = KMC.PRESCRIBED_FORCE
        else:
            index -= 3
            variable = KMC.PRESCRIBED_MOMENT
        values = self.solver.rigid_body_model_part.Nodes[1].GetSolutionStepValue(variable, 0)
        values[index] += prescribed_value
        self.solver.rigid_body_model_part.Nodes[1].SetSolutionStepValue(variable, 0, values)
        #self.solver.prescribed_load[index] += prescribed_value