
import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.rigid_body_process import _RigidBodyProcess


def Factory(solver, parameters):
    return _PrescribedValueProcess(solver, parameters)


class _PrescribedValueProcess(_RigidBodyProcess):

    def __init__(self, solver, parameters):
        super().__init__(solver, parameters)

        self.dof = self.parameters["dof"].GetString()
        self.formula = self.parameters["value"].GetString()
        self.interval = self._AdaptTimeInterval(parameters["interval"])

        if self.dof not in self.solver.available_dofs:
            msg = 'The degree of freedom "' + self.dof + '" is not among the available ones. '
            msg += 'Chose one of the following: ' + str(self.solver.available_dofs)[1:-1] + '.'
            raise Exception(msg)
        elif self.dof not in self.solver.active_dofs:
            msg = 'The degree of freedom "' + self.dof + '" is not active. It needs to'
            msg += 'be activated in the field "active_dofs" in the project parameters.'
            raise Exception(msg)

    def ExecuteInitializeSolutionStep(self):
        if self.solver.time >= self.interval[0] and self.solver.time < self.interval[1]:
            index = self.solver.available_dofs.index(self.dof)
            scope_vars = {'t' : self.solver.time}
            prescribed_value = GenericCallFunction(self.formula, scope_vars, check=False)
            self._ApplyValue(index, prescribed_value)

    def _ApplyValue(self, index, prescribed_value):
        raise Exception("This method should be overwritten in the derived class.")

    def GetDefaultParameters(self):
        return KM.Parameters('''{
                        "dof"      : "UNDEFINED",
                        "interval" : [0, "End"],
                        "value"    : "0"
                    }''')