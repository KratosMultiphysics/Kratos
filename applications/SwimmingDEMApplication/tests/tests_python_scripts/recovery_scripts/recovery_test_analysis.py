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
        self.operator_names = self.scalar_operator_names + self.vector_operator_names

    @staticmethod
    def GetVariableByName(var_name):
        if hasattr(Kratos, var_name):
            return getattr(Kratos, var_name)
        elif hasattr(SDEM, var_name):
            return getattr(SDEM, var_name)
        else:
            raise Exception('Variable '+ var_name + 'has not been adequately added.')


    def SetOperators(self):
        self.scalar_operator_names = []
        self.vector_operator_names = ['scalar_gradient', 'material_derivative']
        GetVariable = RecoveryTestAnalysis.GetVariableByName
        self.vars_man.fluid_vars += [GetVariable('PRESSURE_GRADIENT_ERROR'),
                                     GetVariable('MATERIAL_ACCELERATION_ERROR'),
                                     GetVariable('VELOCITY_DIVERGENCE_ERROR'),
                                     GetVariable('VELOCITY_DIVERGENCE'),
                                     GetVariable('VORTICITY_ERROR'),
                                     GetVariable('VORTICITY')]

    def GetEmbeddedCounter(self):
        return SDP.Counter(is_dead=True)

    def GetDebugInfo(self):
        return SDP.Counter(self.project_parameters["debug_tool_cycle"].GetInt(), 1, is_dead = 1)

    def PerformZeroStepInitializations(self):
        self.local_exact_values = dict()
        self.local_calculated_values = dict()
        self.average_modules = dict()
        self.average_errors = dict()
        self.max_errors = dict()

        for name in self.scalar_operator_names:
            self.local_exact_values[name] = 0.0
            self.local_calculated_values[name] = 0.0
            self.max_errors[name] = 0.0

        for name in self.vector_operator_names:
            self.local_exact_values[name] = Vector([0.0]*3)
            self.local_calculated_values[name] = Vector([0.0]*3)
            self.average_errors[name] = Vector([0.0]*3)

        for name in self.operator_names:
            self.average_modules[name] = 0.0
            self.average_errors[name] = 0.0
            self.max_errors[name] = 0.0

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

    def SetFieldsToImpose(self):
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

    @staticmethod
    def Module(calculated_value):
        if isinstance(calculated_value, float) or isinstance(calculated_value, int):
            return abs(calculated_value)
        elif isinstance(calculated_value, Vector):
            return calculated_value.norm_2()
        else:
            return sum(v**2 for v in list(v)) ** 0.5

    def CalculateRecoveryErrors(self, time):
        total_volume = 0.

        coor = Vector(3)
        error = Vector(3)
        GetVariable = RecoveryTestAnalysis.GetVariableByName

        for node in self.fluid_model_part.Nodes:
            nodal_volume = node.GetSolutionStepValue(Kratos.NODAL_AREA)
            total_volume += nodal_volume

            for i, v in enumerate((node.X, node.Y, node.Z)):
                coor[i] = v

            for name in self.operator_names:

                if name == 'scalar_gradient':
                    self.pressure_field.CalculateGradient(self.time, coor, self.local_exact_values[name], 0)
                    self.local_calculated_values[name] = node.GetSolutionStepValue(GetVariable('PRESSURE_GRADIENT'))
                    error = self.local_calculated_values[name] - self.local_exact_values[name]
                    node.SetSolutionStepValue(GetVariable('PRESSURE_GRADIENT_ERROR'), error)
                elif name == 'material_derivative':
                    self.flow_field.CalculateMaterialAcceleration(self.time, coor, self.local_exact_values[name], 0)
                    self.local_calculated_values[name] = node.GetSolutionStepValue(GetVariable('MATERIAL_ACCELERATION'))
                    error = self.local_calculated_values[name] - self.local_exact_values[name]
                    node.SetSolutionStepValue(GetVariable('MATERIAL_ACCELERATION_ERROR'), error)
                elif name == 'laplacian':
                    self.flow_field.CalculateLaplacian(0., coor, self.local_exact_values[name], 0)
                    self.local_calculated_values[name] = node.GetSolutionStepValue(GetVariable('VELOCITY_LAPLACIAN'))
                    error = self.local_calculated_values[name] - self.local_exact_values[name]
                    node.SetSolutionStepValue(GetVariable('VELOCITY_LAPLACIAN_ERROR'), error)
                elif name == 'divergence':
                    self.local_exact_values[name] = self.flow_field.CalculateDivergence(self.time, coor, 0)
                    self.local_calculated_values[name] = node.GetSolutionStepValue(GetVariable('VELOCITY_DIVERGENCE'))
                    error = self.local_calculated_values[name] - self.local_exact_values[name]
                    node.SetSolutionStepValue(GetVariable('VELOCITY_DIVERGENCE_ERROR'), error)
                elif name == 'rotational':
                    self.flow_field.CalculateRotational(self.time, coor, self.local_exact_values[name], 0)
                    self.local_calculated_values[name] = node.GetSolutionStepValue(GetVariable('VORTICITY'))
                    error = self.local_calculated_values[name] - self.local_exact_values[name]

                    node.SetSolutionStepValue(GetVariable('VORTICITY_ERROR'), error)

                error_norm = RecoveryTestAnalysis.Module(error)
                self.average_modules[name] += RecoveryTestAnalysis.Module(self.local_exact_values[name])  * nodal_volume
                self.average_errors[name] += error_norm  * nodal_volume
                self.max_errors[name] = max(self.max_errors[name], error_norm)

        for name in self.operator_names:
            self.average_errors[name] /= total_volume
            self.average_modules[name] /= total_volume

    def Finalize(self):
        text_width = 30
        super(RecoveryTestAnalysis, self).Finalize()
        text_summary = ''
        text_summary += '\n' + '--' * text_width + '\n'
        for name in self.operator_names:
            text_summary += str(name + ': average modulus').ljust(text_width) + str(self.average_errors[name]) + '\n'
            text_summary += str(name + ': average error').ljust(text_width) + str(self.average_errors[name]) + '\n'
            text_summary += str(name + ': max. error').ljust(text_width) + str(self.max_errors[name]) + '\n'
            text_summary += '--' * text_width + '\n'
        Say(text_summary)