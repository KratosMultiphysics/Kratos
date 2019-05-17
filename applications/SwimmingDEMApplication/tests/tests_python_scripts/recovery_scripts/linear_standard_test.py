import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters, Vector
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from swimming_DEM_analysis import SwimmingDEMAnalysis
from swimming_DEM_analysis import Say
from tests_python_scripts.recovery_scripts.recovery_test_analysis import RecoveryTestAnalysis

class LinearStandardTestAnalysis(RecoveryTestAnalysis):
    def __init__(self, model, varying_parameters=Parameters("{}")):
        super(LinearStandardTestAnalysis, self).__init__(model, varying_parameters)

    def SetOperators(self):
        self.scalar_variable_operator_names = ['gradient']
        self.vector_variable_operator_names = ['divergence', 'rotational', 'gradient']
        GetVariable = RecoveryTestAnalysis.GetVariableByName
        self.vars_man.fluid_vars += [GetVariable('PRESSURE_GRADIENT'),
                                     GetVariable('PRESSURE_GRADIENT_ERROR'),
                                     GetVariable('VELOCITY_DIVERGENCE_ERROR'),
                                     GetVariable('VELOCITY_DIVERGENCE'),
                                     GetVariable('VORTICITY_ERROR'),
                                     GetVariable('SCALAR_GRADIENT'),
                                     GetVariable('VECTOR_GRADIENT'),
                                     GetVariable('VECTOR_GRADIENT_ERROR')]

    def SetFieldsToImpose(self):

        b = 1.0
        bx, by, bz = 1.0, 2.0, 5.0
        b = SDEM.LinearFunction(15.5, b)
        a0 = SDEM.LinearFunction(0.0, bx)
        a1 = SDEM.LinearFunction(0.0, by)
        a2 = SDEM.LinearFunction(0.0, bz)
        self.pressure_field = SDEM.LinearRealField(a0, a1, a2, b)
        field_parameters = Parameters(
            """{
                "A" : [[0.5,1,1], [0,5,-8], [-1,-2,-3]],
                "b" : [1, 0, 0]
                }""")

        self.flow_field = SDEM.LinearVectorField(field_parameters)
        space_time_set = SDEM.SpaceTimeSet()
