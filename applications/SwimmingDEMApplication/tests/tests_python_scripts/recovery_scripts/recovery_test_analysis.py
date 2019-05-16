import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters, Vector, Matrix
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

# According to the Kraots convention:
def GetMatrixIndicesAndKeys():
    indices = []
    keys = []
    for i, i_name in enumerate(('X', 'Y', 'Z')):
        for j, j_name in enumerate(('X', 'Y', 'Z')):
            indices.append((i, j))
            keys.append('_' + i_name + j_name)
    return zip(indices, keys)


class Operator:
    def __init__(self, name, variable_type):
        self.name = name
        self.type = variable_type
        self.derivative_type = self.SetDerivativeType()
        self.InitializeVectors()

    def CalculateNorm(self, variable_type, variable):
        if variable_type == 'scalar':
            return abs(variable)
        elif variable_type == 'vector':
            return variable.norm_2()
        elif variable_type == 'matrix':
            norm = 0.0
            for i in range(3):
                for j in range(3):
                    norm += variable[(i, j)] ** 2

            return norm ** 0.5

    def CalculateErrorNorm(self):
        return float(self.CalculateNorm(self.derivative_type, self.GetError()))

    def CalculateExactValueNorm(self):
        return float(self.CalculateNorm(self.derivative_type, self.local_exact_value))

    def SetDerivativeType(self):
        if self.type == 'scalar':
            if self.name in {'material_derivative', 'laplacian'}:
                return 'scalar'
            elif self.name in {'gradient'}:
                return 'vector'
            else:
                raise Exception('The variable type ' + self.name + 'is not available for variables of type ' + self.type + '.')
        else:
            if self.name in {'divergence'}:
                return 'scalar'
            elif self.name in {'rotational', 'material_derivative', 'laplacian'}:
                return 'vector'
            elif self.name in {'gradient'}:
                return 'matrix'
            else:
                raise Exception('The variable type ' + self.name + 'is not available for variables of type ' + self.type + '.')

    def SetScalarContainers(self):
        self.local_exact_value = 0.0
        self.local_calculated_value = 0.0

    def SetVectorContainers(self):
        self.local_exact_value = Vector([0.0]*3)
        self.local_calculated_value = Vector([0.0]*3)

    def SetMatrixContainers(self):
        self.local_exact_value = Matrix(3,3)
        self.local_calculated_value = Matrix(3,3)

    def SetLocalMatrixValue(self, matrix_value):
        self.local_calculated_value = Matrix(3,3)

        # According to the Kratos convention:
        for i in range(3):
            for j in range(3):
                self.local_calculated_value[(i, j)] = matrix_value[i + j]

    def GetError(self):
        if self.derivative_type in {'scalar', 'vector'}:
            error = self.local_calculated_value - self.local_exact_value
        else:
            error = Matrix(3,3)
            for i in range(3):
                for j in range(3):
                    error[(i,j)] = self.local_calculated_value[(i,j)] - self.local_exact_value[(i,j)]
        return error

    def InitializeVectors(self):
        self.max_error = 0.0
        self.average_module = 0.0
        self.average_error = 0.0

        if self.derivative_type == 'scalar':
            self.SetScalarContainers()
        elif self.derivative_type == 'vector':
            self.SetVectorContainers()
        else:
            self.SetMatrixContainers()

class RecoveryTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters=Parameters("{}")):
        super(RecoveryTestAnalysis, self).__init__(model, varying_parameters)
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'recovery_tests')
        self.project_parameters["custom_fluid"]["fluid_already_calculated"].SetBool(True)
        self.project_parameters.AddEmptyValue("load_derivatives").SetBool(False)
        self.SetOperators()
        self.operators = []
        for name in self.scalar_variable_operator_names:
            self.operators.append(Operator(name, 'scalar'))

        for name in self.vector_variable_operator_names:
            self.operators.append(Operator(name, 'vector'))

    @staticmethod
    def GetVariableByName(var_name):
        if hasattr(Kratos, var_name):
            return getattr(Kratos, var_name)
        elif hasattr(SDEM, var_name):
            return getattr(SDEM, var_name)
        else:
            raise Exception('Variable '+ var_name + 'has not been adequately added.')

    def SetOperators(self):
        self.scalar_variable_operator_names = []
        self.vector_variable_operator_names = ['gradient', 'material_derivative']
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

            for operator in self.operators:
                name = operator.name

                if name == 'gradient':
                    if operator.type == 'scalar':
                        self.pressure_field.CalculateGradient(self.time, coor, operator.local_exact_value, 0)
                        operator.local_calculated_value = node.GetSolutionStepValue(GetVariable('PRESSURE_GRADIENT'))
                        error = operator.GetError()
                        node.SetSolutionStepValue(GetVariable('PRESSURE_GRADIENT_ERROR'), error)
                    else:
                        self.flow_field.CalculateGradient(self.time, coor, operator.local_exact_value, 0)
                        operator.SetLocalMatrixValue(node.GetSolutionStepValue(GetVariable('VECTOR_GRADIENT')))
                        error = operator.GetError()
                        for index_pair, component_key in GetMatrixIndicesAndKeys():
                            node.SetSolutionStepValue(GetVariable('VECTOR_GRADIENT_ERROR' + component_key), error[index_pair])
                elif name == 'material_derivative':
                    self.flow_field.CalculateMaterialAcceleration(self.time, coor, operator.local_exact_value, 0)
                    operator.local_calculated_value = node.GetSolutionStepValue(GetVariable('MATERIAL_ACCELERATION'))
                    error = operator.GetError()
                    node.SetSolutionStepValue(GetVariable('MATERIAL_ACCELERATION_ERROR'), error)
                elif name == 'laplacian':
                    self.flow_field.CalculateLaplacian(0., coor, operator.local_exact_value, 0)
                    operator.local_calculated_value = node.GetSolutionStepValue(GetVariable('VELOCITY_LAPLACIAN'))
                    error = operator.GetError()
                    node.SetSolutionStepValue(GetVariable('VELOCITY_LAPLACIAN_ERROR'), error)
                elif name == 'divergence':
                    operator.local_exact_value = self.flow_field.CalculateDivergence(self.time, coor, 0)
                    operator.local_calculated_value = node.GetSolutionStepValue(GetVariable('VELOCITY_DIVERGENCE'))
                    error = operator.GetError()
                    node.SetSolutionStepValue(GetVariable('VELOCITY_DIVERGENCE_ERROR'), error)
                elif name == 'rotational':
                    self.flow_field.CalculateRotational(self.time, coor, operator.local_exact_value, 0)
                    operator.local_calculated_value = node.GetSolutionStepValue(GetVariable('VORTICITY'))
                    error = operator.GetError()
                    node.SetSolutionStepValue(GetVariable('VORTICITY_ERROR'), error)

                error_norm = operator.CalculateErrorNorm()
                operator.average_module += operator.CalculateExactValueNorm()  * nodal_volume
                operator.average_error += error_norm * nodal_volume
                operator.max_error = max(operator.max_error, error_norm)

        for name in self.operators:
            operator.average_error /= total_volume
            operator.average_module /= total_volume

    def Finalize(self):
        text_width = 30
        super(RecoveryTestAnalysis, self).Finalize()
        text_summary = ''
        text_summary += '\n' + '--' * text_width + '\n'
        for operator in self.operators:
            text_summary += str(operator.name + ': average modulus').ljust(text_width) + str(operator.average_error) + '\n'
            text_summary += str(operator.name + ': average error').ljust(text_width) + str(operator.average_error) + '\n'
            text_summary += str(operator.name + ': max. error').ljust(text_width) + str(operator.max_error) + '\n'
            text_summary += '--' * text_width + '\n'
        Say(text_summary)