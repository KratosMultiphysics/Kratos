import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters, Vector
import swimming_DEM_procedures as SDP
import math
import numpy as np
import os
import sys
import weakref
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from swimming_DEM_analysis import SwimmingDEMAnalysis
from swimming_DEM_analysis import Say
import parameters_tools as PT
import numpy as np

class RecoveryTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters=Parameters("{}")):
        super(RecoveryTestAnalysis, self).__init__(model, varying_parameters)
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'recovery_tests')
        self.project_parameters["custom_fluid"]["fluid_already_calculated"].SetBool(True)
        self.project_parameters.AddEmptyValue("load_derivatives").SetBool(False)
        self.SetOperators()

    def SetOperators(self):
        self.scalar_operator_names = []
        self.vector_operator_names = ['gradient', 'material_derivative']
        self.vars_man.fluid_vars += [SDEM.PRESSURE_GRADIENT_ERROR, SDEM.MATERIAL_ACCELERATION_ERROR]

    def GetEmbeddedCounter(self):
        return SDP.Counter(is_dead=True)

    def GetDebugInfo(self):
        return SDP.Counter(self.project_parameters["debug_tool_cycle"].GetInt(), 1, is_dead = 1)

    def PerformZeroStepInitializations(self):
        self.local_scalar_exact_values = dict()
        self.local_scalar_calculated_values = dict()
        self.scalar_average_modules = dict()
        self.scalar_average_errors = dict()
        self.scalar_max_errors = dict()

        for name in self.scalar_operator_names:
            self.scalar_average_modules[name] = 0.0
            self.scalar_average_errors[name] = 0.0
            self.scalar_max_errors[name] = 0.0

        self.local_vector_exact_values = dict()
        self.local_vector_calculated_values = dict()
        self.vector_average_modules = dict()
        self.vector_average_errors = dict()
        self.vector_max_errors = dict()

        for name in self.vector_operator_names:
            self.local_vector_exact_values[name] = Vector([0.0]*3)
            self.local_vector_calculated_values[name] = Vector([0.0]*3)
            self.vector_average_modules[name] = Vector([0.0]*3)
            self.vector_average_errors[name] = 0.0
            self.vector_max_errors[name] = 0.0

        self.ImposeField()

    def ImposeField(self):
        for node in self.fluid_model_part.Nodes:
            vel= Vector(3)
            coor = Vector([node.X, node.Y, node.Z])
            self.flow_field.Evaluate(self.time, coor, vel, 0)
            pressure = self.pressure_field.Evaluate(self.time, coor, 0)
            node.SetSolutionStepValue(Kratos.VELOCITY, vel)
            node.SetSolutionStepValue(Kratos.PRESSURE, pressure)

    def InitializeSolutionStep(self):
        super(RecoveryTestAnalysis, self).InitializeSolutionStep()
        self.ImposeField()

    def FinalizeSolutionStep(self):
        super(RecoveryTestAnalysis, self).FinalizeSolutionStep()
        self.CalculateRecoveryErrors(self.time)

    def GetFieldUtility(self):
        a = math.pi / 4
        d = math.pi / 2
        b = 1.0
        bx, by, bz = 1.0, 2.0, 5.0
        b = SDEM.LinearFunction(15.5, b)
        a0 = SDEM.LinearFunction(0.0, bx)
        a1 = SDEM.LinearFunction(0.0, by)
        a2 = SDEM.LinearFunction(0.0, bz)
        self.pressure_field = SDEM.LinearRealField(a0, a1, a2, b)
        self.flow_field = SDEM.EthierVelocityField(a, d)
        space_time_set = SDEM.SpaceTimeSet()
        self.field_utility = SDEM.FluidFieldUtility(space_time_set,
                                                    self.pressure_field,
                                                    self.flow_field,
                                                    1000.0,
                                                    1e-6)
        return self.field_utility

    def GetRecoveryCounter(self):
        return SDP.Counter(1, 1, self.project_parameters["coupling"]["coupling_level_type"].GetInt() or self.project_parameters.print_PRESSURE_GRADIENT_option)

    def CalculateRecoveryErrors(self, time):
        total_volume = 0.

        coor = Vector(3)
        error = Vector(3)

        for node in self.fluid_model_part.Nodes:
            nodal_volume = node.GetSolutionStepValue(Kratos.NODAL_AREA)
            total_volume += nodal_volume

            for i, v in enumerate((node.X, node.Y, node.Z)):
                coor[i] = v

            for name in self.vector_operator_names:

                if name == 'gradient':
                    self.pressure_field.CalculateGradient(self.time, coor, self.local_vector_exact_values[name], 0)
                    self.local_vector_calculated_values[name] = node.GetSolutionStepValue(Kratos.PRESSURE_GRADIENT)
                    error = self.local_vector_calculated_values[name] - self.local_vector_exact_values[name]
                    node.SetSolutionStepValue(SDEM.PRESSURE_GRADIENT_ERROR, error)
                elif name == 'material_derivative':
                    self.flow_field.CalculateMaterialAcceleration(self.time, coor, self.local_vector_exact_values[name], 0)
                    self.local_vector_calculated_values[name] = node.GetSolutionStepValue(Kratos.MATERIAL_ACCELERATION)
                    error = self.local_vector_calculated_values[name] - self.local_vector_exact_values[name]
                    node.SetSolutionStepValue(SDEM.MATERIAL_ACCELERATION_ERROR, error)
                elif name == 'laplacian':
                    self.flow_field.CalculateLaplacian(0., coor, self.local_vector_exact_values[name], 0)
                    self.local_vector_calculated_values[name] = node.GetSolutionStepValue(Kratos.VELOCITY_LAPLACIAN)
                    error = self.local_vector_calculated_values[name] - self.local_vector_exact_values[name]
                    node.SetSolutionStepValue(SDEM.VELOCITY_LAPLACIAN_ERROR, error)

                error_norm = error.norm_2()
                self.vector_average_modules[name] += self.local_vector_exact_values[name].norm_2()  * nodal_volume
                self.vector_average_errors[name] += error_norm  * nodal_volume
                self.vector_max_errors[name] = max(self.vector_max_errors[name], error_norm)

        self.vector_average_errors[name] /= total_volume
        self.vector_average_modules[name] /= total_volume

    def Finalize(self):
        text_width = 45
        super(RecoveryTestAnalysis, self).Finalize()
        Say('\n' + '-.' * text_width)
        for name in self.vector_operator_names:
            if self.vector_average_errors[name] > 0:
                Say(str('Average error for the ' + name).ljust(text_width), self.vector_average_errors[name])
                Say(str('Max. error for the ' + name).ljust(text_width), self.vector_max_errors[name])
        Say('-.' * text_width + '\n')