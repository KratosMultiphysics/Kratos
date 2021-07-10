import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from math import isclose


class TestResponseFunction(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

        prop = cls.model_part.GetProperties()[0]
        prop[Kratos.DENSITY] = 1.5
        prop[Kratos.DYNAMIC_VISCOSITY] = 1.2e-4
        cls.model_part.CreateNewElement("Element2D3N", 1, [2, 3, 1], prop)

        cls.model_part.SetBufferSize(1)

        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VELOCITY, 1.0, 100.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.ACCELERATION, 0.0, 50.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.PRESSURE, 0.0, 10.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.BODY_FORCE, 0.0, 20.0, 0)

        cls.model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.04
        cls.model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 0.3

        cls.response_function = KratosCFD.ResidualResponseFunction2D(Kratos.Parameters("""{
            "continuity_residual_weight": 15.0,
            "momentum_residual_weight": 24.0
        }"""),  cls.model_part)
        cls.ref_value = cls.response_function.CalculateValue(cls.model_part)
        cls.element = cls.model_part.GetElement(1)
        cls.domain_size = 2
        cls.number_of_nodes = 3
        cls.block_size = 3

    def testCalculateFirstDerivativesGradient(self):
        residual_local_size = self.block_size * self.number_of_nodes
        delta = 1e-5

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        self.response_function.CalculateFirstDerivativesGradient(
            self.element, residual_gradient, response_gradient, self.model_part.ProcessInfo)

        fd_response_gradient = Kratos.Vector(response_gradient.Size(), 0.0)

        for c, node in enumerate(self.element.GetGeometry()):
            for k in range(self.domain_size):
                velocity = node.GetSolutionStepValue(Kratos.VELOCITY)

                velocity[k] += delta
                node.SetSolutionStepValue(Kratos.VELOCITY, velocity)

                fd_response_gradient[c * self.block_size + k] = self._CalculateSensitivity(delta)

                velocity[k] -= delta
                node.SetSolutionStepValue(Kratos.VELOCITY, velocity)

            pressure = node.GetSolutionStepValue(Kratos.PRESSURE)

            pressure += delta
            node.SetSolutionStepValue(Kratos.PRESSURE, pressure)

            fd_response_gradient[c * self.block_size + self.domain_size] = self._CalculateSensitivity(delta)

            pressure -= delta
            node.SetSolutionStepValue(Kratos.PRESSURE, pressure)

        self._IsVectorRelativelyClose(fd_response_gradient, response_gradient, 1e-5, 1e-5)

    def testCalculateSecondDerivativesGradient(self):
        residual_local_size = self.block_size * self.number_of_nodes
        delta = 2e-5

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        self.response_function.CalculateSecondDerivativesGradient(
            self.element, residual_gradient, response_gradient, self.model_part.ProcessInfo)

        fd_response_gradient = Kratos.Vector(response_gradient.Size(), 0.0)

        for c, node in enumerate(self.element.GetGeometry()):
            for k in range(self.domain_size):
                acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)

                acceleration[k] += delta
                node.SetSolutionStepValue(Kratos.ACCELERATION, acceleration)

                fd_response_gradient[c * self.block_size + k] = self._CalculateSensitivity(delta)

                acceleration[k] -= delta
                node.SetSolutionStepValue(Kratos.ACCELERATION, acceleration)

        self._IsVectorRelativelyClose(fd_response_gradient, response_gradient, 1e-6, 1e-5)

    def testCalculatePartialSensitivity(self):
        residual_local_size = self.domain_size * self.number_of_nodes
        delta = 1e-7

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        self.response_function.CalculatePartialSensitivity(
            self.element, Kratos.SHAPE_SENSITIVITY, residual_gradient, response_gradient, self.model_part.ProcessInfo)

        fd_response_gradient = Kratos.Vector(response_gradient.Size(), 0.0)

        for c, node in enumerate(self.element.GetGeometry()):
            node.X += delta
            fd_response_gradient[c * self.domain_size + 0] = self._CalculateSensitivity(delta)
            node.X -= delta

            node.Y += delta
            fd_response_gradient[c * self.domain_size + 1] = self._CalculateSensitivity(delta)
            node.Y -= delta

        self._IsVectorRelativelyClose(fd_response_gradient, response_gradient, 1e-6, 1e-5)

    def _CalculateSensitivity(self, delta):
        value = self.response_function.CalculateValue(self.model_part)
        return (value - self.ref_value) / delta

    def _IsVectorRelativelyClose(self, vec_a, vec_b, rel_tol, abs_tol):
        self.assertEqual(vec_a.Size(), vec_b.Size())
        for i in range(vec_a.Size()):
            if (not isclose(vec_b[i], vec_a[i], rel_tol=rel_tol, abs_tol=abs_tol)):
                msg = "VecA[{:d}] != VecB[{:d}] [ {:f} != {:f} ]. Vectors are: \n\t VecA = {:s}\n\t VecB = {:s}".format(
                    i, i, vec_a[i], vec_b[i], str(vec_a), str(vec_b))
                raise AssertionError(msg)

class TestDomainIntegratedResponseFunction(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

        prop = cls.model_part.GetProperties()[0]
        prop[Kratos.DENSITY] = 1.5
        prop[Kratos.DYNAMIC_VISCOSITY] = 1.2e-4
        cls.model_part.CreateNewElement("Element2D3N", 1, [2, 3, 1], prop)

        cls.model_part.SetBufferSize(1)

        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.PRESSURE, 1.0, 100.0, 0)

        cls.model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.04
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        cls.response_function = KratosCFD.DomainIntegratedResponseFunction(Kratos.Parameters("""{
            "variable_name": "PRESSURE"
        }"""),  cls.model_part)

        Kratos.VariableUtils().SetFlag(Kratos.STRUCTURE, True, cls.model_part.Elements)
        cls.ref_value = cls.response_function.CalculateValue(cls.model_part)
        cls.element = cls.model_part.GetElement(1)
        cls.domain_size = 2
        cls.number_of_nodes = 3
        cls.block_size = 1

    def testCalculateFirstDerivativesGradient(self):
        residual_local_size = self.block_size * self.number_of_nodes
        delta = 1e-5

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        self.response_function.CalculateFirstDerivativesGradient(
            self.element, residual_gradient, response_gradient, self.model_part.ProcessInfo)

        fd_response_gradient = Kratos.Vector(response_gradient.Size(), 0.0)

        for c, node in enumerate(self.element.GetGeometry()):
            pressure = node.GetSolutionStepValue(Kratos.PRESSURE)

            pressure += delta
            node.SetSolutionStepValue(Kratos.PRESSURE, pressure)

            fd_response_gradient[c * self.block_size] = self._CalculateSensitivity(delta)

            pressure -= delta
            node.SetSolutionStepValue(Kratos.PRESSURE, pressure)

        self._IsVectorRelativelyClose(fd_response_gradient, response_gradient, 1e-5, 1e-5)

    def testCalculateSecondDerivativesGradient(self):
        residual_local_size = self.block_size * self.number_of_nodes
        delta = 2e-5

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        self.response_function.CalculateSecondDerivativesGradient(
            self.element, residual_gradient, response_gradient, self.model_part.ProcessInfo)

        fd_response_gradient = Kratos.Vector(response_gradient.Size(), 0.0)

        self._IsVectorRelativelyClose(fd_response_gradient, response_gradient, 1e-6, 1e-5)

    def testCalculatePartialSensitivity(self):
        residual_local_size = self.domain_size * self.number_of_nodes
        delta = 1e-7

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        self.response_function.CalculatePartialSensitivity(
            self.element, Kratos.SHAPE_SENSITIVITY, residual_gradient, response_gradient, self.model_part.ProcessInfo)

        fd_response_gradient = Kratos.Vector(response_gradient.Size(), 0.0)

        for c, node in enumerate(self.element.GetGeometry()):
            node.X += delta
            fd_response_gradient[c * self.domain_size + 0] = self._CalculateSensitivity(delta)
            node.X -= delta

            node.Y += delta
            fd_response_gradient[c * self.domain_size + 1] = self._CalculateSensitivity(delta)
            node.Y -= delta

        self._IsVectorRelativelyClose(fd_response_gradient, response_gradient, 1e-6, 1e-5)

    def _CalculateSensitivity(self, delta):
        value = self.response_function.CalculateValue(self.model_part)
        return (value - self.ref_value) / delta

    def _IsVectorRelativelyClose(self, vec_a, vec_b, rel_tol, abs_tol):
        self.assertEqual(vec_a.Size(), vec_b.Size())
        for i in range(vec_a.Size()):
            if (not isclose(vec_b[i], vec_a[i], rel_tol=rel_tol, abs_tol=abs_tol)):
                msg = "VecA[{:d}] != VecB[{:d}] [ {:f} != {:f} ]. Vectors are: \n\t VecA = {:s}\n\t VecB = {:s}".format(
                    i, i, vec_a[i], vec_b[i], str(vec_a), str(vec_b))
                raise AssertionError(msg)

if __name__ == '__main__':
    UnitTest.main()
