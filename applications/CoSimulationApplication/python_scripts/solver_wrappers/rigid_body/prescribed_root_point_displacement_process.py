
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.prescribed_value_process import _PrescribedValueProcess


def Factory(solver, parameters):
    return PrescribedRootPointDisplacementProcess(solver, parameters)


class PrescribedRootPointDisplacementProcess(_PrescribedValueProcess):
    
    def _ApplyValue(self, index, prescribed_value):
        self.solver.prescribed_root_point_displ[index] += prescribed_value