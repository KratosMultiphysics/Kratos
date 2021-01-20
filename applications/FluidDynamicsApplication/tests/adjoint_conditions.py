import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.KratosUnittest as KratosUnittest

import random

class TestAdjointMonolithicWallCondition(KratosUnittest.TestCase):
    def setUp(self):
        self.delta_time = 1.0
        # create test model part
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)

        self.model_part.ProcessInfo[Kratos.OSS_SWITCH] = 0
        self.model_part.ProcessInfo[Kratos.DELTA_TIME] = self.delta_time
        self.model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 1.0

        self.model_part.SetBufferSize(2)

        prop = self.model_part.GetProperties()[0]
        self.model_part.CreateNewCondition("MonolithicWallCondition2D2N", 1, [1, 2], prop)
        self.model_part.CreateNewCondition("AdjointMonolithicWallCondition2D2N", 2, [1, 2], prop)

        self.primal_condition = self.model_part.GetCondition(1)
        self.adjoint_condition = self.model_part.GetCondition(2)

        self._AssignSolutionStepData(0)

        Kratos.NormalCalculationUtils().CalculateOnSimplex(self.model_part, 2)
        Kratos.NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(self.model_part.Conditions, 2)

    def _AssignSolutionStepData(self, step=0):
        # generate nodal solution step test data
        random.seed(1.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DENSITY, step, 1.0)
            node.SetSolutionStepValue(Kratos.VISCOSITY, step, 1.0e-5)
            node.SetSolutionStepValue(Kratos.VELOCITY_X, step, random.random())
            node.SetSolutionStepValue(Kratos.VELOCITY_Y, step, random.random())
            node.SetSolutionStepValue(Kratos.ACCELERATION_X, step, random.random())
            node.SetSolutionStepValue(Kratos.ACCELERATION_Y, step, random.random())
            node.SetSolutionStepValue(Kratos.PRESSURE, step, random.random())
            node.SetSolutionStepValue(Kratos.EXTERNAL_PRESSURE, step, random.random())
            node.SetValue(KratosCFD.Y_WALL, random.random())

    def _SetSlip(self):
        Kratos.VariableUtils().SetFlag(Kratos.SLIP, True, self.model_part.Nodes)
        Kratos.VariableUtils().SetFlag(Kratos.SLIP, True, self.model_part.Conditions)

    def _assertVectorAlmostEqual(self, vector1, vector2, prec=7):
        self.assertEqual(vector1.Size(), vector2.Size())
        for i in range(vector1.Size()):
            self.assertAlmostEqual(vector1[i], vector2[i], prec)

    def _assertMatrixAlmostEqual(self, matrix1, matrix2, prec=7):
        self.assertEqual(matrix1.Size1(), matrix2.Size1())
        self.assertEqual(matrix1.Size2(), matrix2.Size2())
        for i in range(matrix1.Size1()):
            for j in range(matrix1.Size2()):
                self.assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

    def _CalculateResidual(self, condition):
        lhs = Kratos.Matrix()
        rhs = Kratos.Vector()
        condition.CalculateLocalSystem(lhs, rhs, self.model_part.ProcessInfo)
        condition.CalculateLocalVelocityContribution(lhs, rhs, self.model_part.ProcessInfo)

        return rhs

    def _GetPerturbationMethod(self, variable):
        if (variable == Kratos.SHAPE_SENSITIVITY):
            def perturbation_method(node, index, value):
                if (index == 0):
                    node.X += value
                elif (index == 1):
                    node.Y += value
                else:
                    node.Z += value
            self.dimensionality = 2
        elif (Kratos.KratosGlobals.GetVariableType(variable.Name()) == "Array"):
            def perturbation_method(node, index, value):
                var_value = node.GetSolutionStepValue(variable)
                var_value[index] += value
                node.SetSolutionStepValue(variable, 0, var_value)
            self.dimensionality = 2
        elif (Kratos.KratosGlobals.GetVariableType(variable.Name()) == "Double"):
            def perturbation_method(node, index, value):
                var_value = node.GetSolutionStepValue(variable)
                var_value += value
                node.SetSolutionStepValue(variable, 0, var_value)
            self.dimensionality = 1
        else:
            def perturbation_method(node, index, value):
                raise Exception("Unknown variable type")

        return perturbation_method

    def _ComputeFiniteDifferenceSensitivity(self, variables_list, delta):
        residual_0 = self._CalculateResidual(self.primal_condition)
        total_dimensionality = 0
        for variable in variables_list:
            _ = self._GetPerturbationMethod(variable)
            total_dimensionality += self.dimensionality
        fd_sensitivity = Kratos.Matrix(len(self.primal_condition.GetGeometry()) * total_dimensionality, residual_0.Size())
        local_index = 0
        for node in self.primal_condition.GetGeometry():
            for variable in variables_list:
                perturbation_method = self._GetPerturbationMethod(variable)
                for i in range(self.dimensionality):
                    perturbation_method(node, i, delta)

                    residual = self._CalculateResidual(self.primal_condition)
                    for j in range(residual.Size()):
                        fd_sensitivity[local_index, j] = (residual[j] - residual_0[j]) / delta

                    local_index += 1

                    perturbation_method(node, i, -delta)

        return fd_sensitivity

    def testResiduals(self):
        primal_residual = self._CalculateResidual(self.primal_condition)
        adjoint_residual = self._CalculateResidual(self.adjoint_condition)

        self._assertVectorAlmostEqual(primal_residual, adjoint_residual, 9)

    def testResidualsSlip(self):
        self._SetSlip()
        primal_residual = self._CalculateResidual(self.primal_condition)
        adjoint_residual = self._CalculateResidual(self.adjoint_condition)

        self._assertVectorAlmostEqual(primal_residual, adjoint_residual, 9)

    def testCalculateSensitivityMatrix(self):
        analytical_sensitivity = Kratos.Matrix()
        self.adjoint_condition.CalculateSensitivityMatrix(Kratos.SHAPE_SENSITIVITY, analytical_sensitivity, self.model_part.ProcessInfo)
        fd_sensitivity = self._ComputeFiniteDifferenceSensitivity([Kratos.SHAPE_SENSITIVITY], 1e-6)
        self._assertMatrixAlmostEqual(fd_sensitivity, analytical_sensitivity, 9)

    def testCalculateSensitivityMatrixSlip(self):
        self._SetSlip()
        analytical_sensitivity = Kratos.Matrix()
        self.adjoint_condition.CalculateSensitivityMatrix(Kratos.SHAPE_SENSITIVITY, analytical_sensitivity, self.model_part.ProcessInfo)
        fd_sensitivity = self._ComputeFiniteDifferenceSensitivity([Kratos.SHAPE_SENSITIVITY], 1e-7)
        self._assertMatrixAlmostEqual(fd_sensitivity, analytical_sensitivity, 9)

    def testCalculateSecondDerivativesLHS(self):
        analytical_sensitivity = Kratos.Matrix()
        self.adjoint_condition.CalculateSecondDerivativesLHS(analytical_sensitivity, self.model_part.ProcessInfo)
        self._assertMatrixAlmostEqual(Kratos.Matrix(6, 6, 0), analytical_sensitivity, 9)

    def testCalculateSecondDerivativesLHSSlip(self):
        self._SetSlip()
        analytical_sensitivity = Kratos.Matrix()
        self.adjoint_condition.CalculateSecondDerivativesLHS(analytical_sensitivity, self.model_part.ProcessInfo)
        self._assertMatrixAlmostEqual(Kratos.Matrix(6, 6, 0), analytical_sensitivity, 9)

    def testCalculateFirstDerivativesLHS(self):
        analytical_sensitivity = Kratos.Matrix()
        self.adjoint_condition.CalculateFirstDerivativesLHS(analytical_sensitivity, self.model_part.ProcessInfo)
        self._assertMatrixAlmostEqual(Kratos.Matrix(6, 6, 0), analytical_sensitivity, 9)

    def testCalculateFirstDerivativesLHSSlip(self):
        self._SetSlip()
        analytical_sensitivity = Kratos.Matrix()
        self.adjoint_condition.CalculateFirstDerivativesLHS(analytical_sensitivity, self.model_part.ProcessInfo)
        fd_sensitivity = self._ComputeFiniteDifferenceSensitivity([Kratos.VELOCITY, Kratos.PRESSURE], 1e-7)
        self._assertMatrixAlmostEqual(fd_sensitivity, analytical_sensitivity, 9)

if __name__ == '__main__':
    KratosUnittest.main()
