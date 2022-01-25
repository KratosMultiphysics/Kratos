
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.prescribed_value_process import _PrescribedValueProcess


def Factory(solver, parameters):
    return PrescribedForceProcess(solver, parameters)


class PrescribedForceProcess(_PrescribedValueProcess):
    
    def _ApplyValue(self, prescribed_value):
        self.solver.prescribed_load += prescribed_value