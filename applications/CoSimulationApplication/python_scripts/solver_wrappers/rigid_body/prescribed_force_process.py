
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.prescribed_value_process import _PrescribedValueProcess


def Factory(solver, parameters):
    return PrescribedForceProcess(solver, parameters)


class PrescribedForceProcess(_PrescribedValueProcess):
    
    def _ApplyValue(self, index, prescribed_value):
        self.solver.prescribed_load[index] += prescribed_value