
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.rigid_body_process import _RigidBodyProcess


def Factory(solver, parameters):
    return _PrescribedValueProcess(solver, parameters)


class SaveRestartProcess(_RigidBodyProcess):

    def __init__(self, solver, parameters):
        super().__init__(solver, parameters)