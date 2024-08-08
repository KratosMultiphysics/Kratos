import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestLSSVariableUtilities(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # create test model part
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        # add solution step variables
        variables_list = [
            Kratos.VELOCITY,
            Kratos.PRESSURE,
            Kratos.ACCELERATION,
            KratosCFD.LSS_VELOCITY,
            KratosCFD.LSS_PRESSURE,
            KratosCFD.LSS_ACCELERATION,
            KratosCFD.ADJOINT_FLUID_VECTOR_1,
            KratosCFD.ADJOINT_FLUID_SCALAR_1,
            KratosCFD.ADJOINT_FLUID_VECTOR_3
        ]

        for var in variables_list:
            cls.model_part.AddNodalSolutionStepVariable(var)

        cls.model_part.SetBufferSize(2)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

        prop = cls.model_part.GetProperties()[0]
        cls.condition = cls.model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], prop)
        cls.element = cls.model_part.CreateNewElement("VMS2D3N", 1, [1, 2, 3], prop)

        for var in variables_list:
            KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, var, 5.0, 100.0, 0)
            KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, var, 10.0, 1000.0, 1)

        cls.primal_variables_list = [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.PRESSURE]
        cls.primal_first_derivative_variables_list = [Kratos.ACCELERATION_X, Kratos.ACCELERATION_Y, None]

        cls.lss_variables_list = [KratosCFD.LSS_VELOCITY_X, KratosCFD.LSS_VELOCITY_Y, KratosCFD.LSS_PRESSURE]
        cls.lss_first_derivative_variables_list = [KratosCFD.LSS_ACCELERATION_X, KratosCFD.LSS_ACCELERATION_Y, None]

        cls.adjoint_variables_list = [KratosCFD.ADJOINT_FLUID_VECTOR_1_X, KratosCFD.ADJOINT_FLUID_VECTOR_1_Y, KratosCFD.ADJOINT_FLUID_SCALAR_1]
        cls.adjoint_first_derivative_variables_list = [KratosCFD.ADJOINT_FLUID_VECTOR_3_X, KratosCFD.ADJOINT_FLUID_VECTOR_3_Y, None]

        cls.fluid_lss_variable_utilities = KratosCFD.FluidLSSVariableUtilities(
            cls.primal_variables_list,
            cls.primal_first_derivative_variables_list,
            cls.adjoint_variables_list,
            cls.adjoint_first_derivative_variables_list,
            cls.lss_variables_list,
            cls.lss_first_derivative_variables_list
        )

    def testGetPrimalValues(self):
        self.__CheckMethod(self.fluid_lss_variable_utilities.GetPrimalValues, self.primal_variables_list)

    def testGetPrimalFirstDerivativeValues(self):
        self.__CheckMethod(self.fluid_lss_variable_utilities.GetPrimalFirstDerivativeValues, self.primal_first_derivative_variables_list)

    def testGetAdjointValues(self):
        self.__CheckMethod(self.fluid_lss_variable_utilities.GetAdjointValues, self.adjoint_variables_list)

    def testGetAdjointFirstDerivativeValues(self):
        self.__CheckMethod(self.fluid_lss_variable_utilities.GetAdjointFirstDerivativeValues, self.adjoint_first_derivative_variables_list)

    def testGetLSSValues(self):
        self.__CheckMethod(self.fluid_lss_variable_utilities.GetLSSValues, self.lss_variables_list)

    def testGetLSSFirstDerivativeValues(self):
        self.__CheckMethod(self.fluid_lss_variable_utilities.GetLSSFirstDerivativeValues, self.lss_first_derivative_variables_list)

    def __CheckMethod(self, method, variables_list):
        current_v = Kratos.Vector()
        method(current_v, self.condition, 0)
        self.__CheckVector(self.condition, variables_list, 0, current_v)

        old_v = Kratos.Vector()
        method(old_v, self.condition, 1)
        self.__CheckVector(self.condition, variables_list, 1, old_v)

        current_v = Kratos.Vector()
        method(current_v, self.element, 0)
        self.__CheckVector(self.element, variables_list, 0, current_v)

        old_v = Kratos.Vector()
        method(old_v, self.element, 1)
        self.__CheckVector(self.element, variables_list, 1, old_v)

    def __CheckVector(self, entity, variables_list, step, values):
        index = 0
        for node in entity.GetGeometry():
            for var in variables_list:
                if var is not None:
                    self.assertAlmostEqual(values[index], node.GetSolutionStepValue(var, step), 12)
                else:
                    self.assertEqual(values[index], 0.0)
                index += 1

if __name__ == '__main__':
    KratosUnittest.main()
