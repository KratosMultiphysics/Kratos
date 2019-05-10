import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters, Vector
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from swimming_DEM_analysis import SwimmingDEMAnalysis
from swimming_DEM_analysis import Say
from tests_python_scripts.recovery_scripts.recovery_test_analysis import RecoveryTestAnalysis

class EthierStandardTestAnalysis(RecoveryTestAnalysis):
    def __init__(self, model, varying_parameters=Parameters("{}")):
        super(EthierStandardTestAnalysis, self).__init__(model, varying_parameters)

    def SetOperators(self):
        self.scalar_operator_names = []
        self.vector_operator_names = ['gradient', 'material_derivative']
        self.vars_man.fluid_vars += [SDEM.PRESSURE_GRADIENT_ERROR, SDEM.MATERIAL_ACCELERATION_ERROR]

    def GetFieldUtility(self):
        import math
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
