'''HDF5 controllers.

This module contains the controllers which control the frequency of HDF5 IO
operations.

license: HDF5Application/license.txt
'''


import KratosMultiphysics


class DefaultController(object):
    '''Simple pass through controller.'''

    def __init__(self, model_part, io, settings=None):
        self.model_part = model_part
        self.io = io
        self.operations = []

    def Add(self, operation):
        self.operations.append(operation)

    def __call__(self):
        current_io = self.io.Get(self.model_part)
        for op in self.operations:
            op(self.model_part, current_io)


class TemporalController(object):
    '''Frequency-based controller.

    Controls execution according to the 'time_frequency' and 'step_frequency'
    specified in the json settings.
    '''

    def __init__(self, model_part, io, settings):
        self.model_part = model_part
        self.io = io
        settings.SetDefault('time_frequency', 1.0)
        settings.SetDefault('step_frequency', 1)
        self.time_frequency = settings['time_frequency']
        self.step_frequency = settings['step_frequency']
        self.operations = []
        self.current_time = 0.0
        self.current_step = 0

    def Add(self, operation):
        self.operations.append(operation)

    def _IsOutputStep(self):
        if self.current_step == self.step_frequency:
            return True
        if self.current_time > self.time_frequency:
            return True
        # Compare relative error against an epsilon, which is much larger than
        # the machine epsilon, and include a lower bound based on
        # https://github.com/chromium/chromium, cc::IsNearlyTheSame.
        eps = 1e-6
        tol = eps * max(abs(self.current_time), abs(self.time_frequency), eps)
        if abs(self.current_time - self.time_frequency) < tol:
            return True
        return False

    def __call__(self):
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.current_time += delta_time
        self.current_step += 1
        if self._IsOutputStep():
            current_io = self.io.Get(self.model_part)
            for op in self.operations:
                op(self.model_part, current_io)
            self.current_time = 0.0
            self.current_step = 0


def Create(model_part, io, settings):
    '''Return the controller specified by the setting 'controller_type'.

    Empty settings will contain default values after returning from the
    function call.
    '''
    settings.SetDefault('controller_type', 'default_controller')
    controller_type = settings['controller_type']
    if controller_type == 'default_controller':
        return DefaultController(model_part, io, settings)
    elif controller_type == 'temporal_controller':
        return TemporalController(model_part, io, settings)
    else:
        raise ValueError(
            '"controller_type" has invalid value "' + controller_type + '"')
