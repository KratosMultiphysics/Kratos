
import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction

import numpy as np

def CreateRigidBodyProcess(solver, settings):

    if settings.Has("name"):
        process_name = settings["name"].GetString()
        parameters = settings["parameters"]

        if process_name == "impose_prescribed_force":
            return PrescribedForceProcess(solver, parameters)

        elif process_name == "impose_prescribed_root_point_displacement":
            return PrescribedRootPointDisplacementProcess(solver, parameters)

        else:
            msg = 'The process "' + process_name + '" for the rigid body solver '
            msg = 'is not among the available ones. Choose a valid boundary condition type.'
            raise Exception(msg)
    else:
        msg = 'The field "name" must be provided when defining a RigidBodyProcess.'
        raise Exception(msg)


class _RigidBodyProcess():

    def __init__(self, solver, parameters):

        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())
        self.parameters = parameters
        self.solver = solver

    def ExecuteInitialize(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass
    
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteFinalizeSolutionStep(self):
        pass
    
    def ExecuteBeforeOutputStep(self):
        pass
    
    def ExecuteAfterOutputStep(self):
        pass
    
    def ExecuteFinalize(self):
        pass

    def GetDefaultParameters(self):
        pass
    
    #TODO: how to check a KM.Parameters and loop through it
    # Maybe check ValidateAndAssign... or assign_scalar_*_process
    def _AdaptTimeInterval(self, interval):
        if interval[1].IsString():
            if interval[1].GetString() == "End":
                interval[1].SetDouble(1e30)
        interval = np.array([interval[0].GetDouble(), interval[1].GetDouble()])
        return interval


class _PrescrivedValueProcess(_RigidBodyProcess):

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

    def ExecuteBeforeSolutionLoop(self):
        if self.solver.time >= self.interval[0] and self.solver.time < self.interval[1]:
            index = self.solver.available_dofs.index(self.dof)
            scope_vars = {'t' : self.solver.time}
            prescribed_value = GenericCallFunction(self.formula, scope_vars, check=False)
            self._ApplyValue(prescribed_value)

    def _ApplyValue(self, prescribed_value):
        pass

    def GetDefaultParameters(self):
        return KM.Parameters('''{
                        "dof"      : "UNDEFINED",
                        "interval" : [0, "End"],
                        "value"    : "0"
                    }''')

class PrescribedForceProcess(_PrescrivedValueProcess):
    
    def _ApplyValue(self, prescribed_value):
        self.solver.prescribed_load += prescribed_value

class PrescribedRootPointDisplacementProcess(_PrescrivedValueProcess):
    
    def _ApplyValue(self, prescribed_value):
        self.solver.prescribed_root_point_displ += prescribed_value


    
    