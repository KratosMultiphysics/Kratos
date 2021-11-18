import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from math import isclose, pow


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


class TestDomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction3D(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.sub_model_part = cls.model_part.CreateSubModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        cls.sub_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.sub_model_part.CreateNewNode(2, 1.0, 0.0, 1.5)
        cls.sub_model_part.CreateNewNode(3, 1.4, 1.0, 0.3)
        cls.sub_model_part.CreateNewNode(4, 0.7, 0.4, 1.8)

        prop = cls.model_part.GetProperties()[0]
        prop[Kratos.DENSITY] = 1.5
        prop[Kratos.DYNAMIC_VISCOSITY] = 1.2e-4
        cls.sub_model_part.CreateNewElement("Element2D3N", 1, [2, 3, 1], prop)
        cls.sub_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [2, 3, 1], prop)
        cls.sub_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [2, 3, 4], prop)

        cls.eq_id_map = {
            1: [1, 2, 0],
            2: [1, 2, 3],
        }

        cls.model_part.SetBufferSize(1)

        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VELOCITY, 1.0, 100.0, 0)

        cls.model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.04
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.value_to_power = 2

        parameters = Kratos.Parameters("""{
            "model_part_name"          : "test",
            "variable_name"            : "VELOCITY",
            "magnitude_square_to_power": 2
        }""")

        cls.response_function = KratosCFD.DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction3D(parameters,  cls.model_part)

        cls.response_function.Initialize()
        cls.response_function.InitializeSolutionStep()
        cls.ref_value = cls.response_function.CalculateValue(cls.model_part)
        cls.condition = cls.model_part.GetCondition(1)
        cls.domain_size = 3
        cls.number_of_nodes = 3
        cls.block_size = 3

    def testCalculateValue(self):
        value = 0.0
        total_area = 0.0
        for condition in self.model_part.Conditions:
            area = condition.GetGeometry().DomainSize()
            velocity = Kratos.Array3(0.0)
            for node in condition.GetGeometry():
                velocity += node.GetSolutionStepValue(Kratos.VELOCITY) / 3.0
            value += area * pow((velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]), self.value_to_power)
            total_area += area

        self.assertAlmostEqual(value / total_area, self.ref_value, 9)

    def testCalculateFirstDerivativesGradient(self):
        residual_local_size = self.block_size * self.number_of_nodes
        delta = 1e-8

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        fd_response_gradient = Kratos.Vector(self.block_size * self.model_part.NumberOfNodes(), 0.0)
        assembled_adjoint_gradient = Kratos.Vector(self.block_size * self.model_part.NumberOfNodes(), 0.0)

        for condition in self.model_part.Conditions:
            self.response_function.CalculateFirstDerivativesGradient(
                condition, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            self._AssembleAdjointGradients(condition, assembled_adjoint_gradient, response_gradient)

        for c, node in enumerate(self.model_part.Nodes):
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            for k in range(3):
                velocity[k] += delta
                node.SetSolutionStepValue(Kratos.VELOCITY, velocity)

                fd_response_gradient[c * self.block_size + k] = self._CalculateSensitivity(delta)

                velocity[k] -= delta
                node.SetSolutionStepValue(Kratos.VELOCITY, velocity)

        self._IsVectorRelativelyClose(fd_response_gradient, assembled_adjoint_gradient, 1e-5, 1e-5)

    def testCalculateSecondDerivativesGradient(self):
        residual_local_size = self.block_size * self.number_of_nodes

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        self.response_function.CalculateSecondDerivativesGradient(
            self.condition, residual_gradient, response_gradient, self.model_part.ProcessInfo)

        self._IsVectorRelativelyClose(Kratos.Vector(response_gradient.Size(), 0.0), response_gradient, 1e-5, 1e-5)

    def testCalculatePartialSensitivity(self):
        residual_local_size = self.block_size * self.number_of_nodes
        delta = 1e-5

        response_gradient = Kratos.Vector()
        residual_gradient = Kratos.Matrix(residual_local_size, residual_local_size)

        fd_response_gradient = Kratos.Vector(self.block_size * self.model_part.NumberOfNodes(), 0.0)
        assembled_adjoint_gradient = Kratos.Vector(self.block_size * self.model_part.NumberOfNodes(), 0.0)

        def node_coordinates_setter(node, direction, delta):
            if direction == 0:
                node.X = node.X + delta
            elif direction == 1:
                node.Y = node.Y + delta
            elif direction == 2:
                node.Z = node.Z + delta
            else:
                raise Exception("Unsupported direction")

        for condition in self.model_part.Conditions:
            self.response_function.CalculatePartialSensitivity(
                condition, Kratos.SHAPE_SENSITIVITY, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            self._AssembleAdjointGradients(condition, assembled_adjoint_gradient, response_gradient)

        for c, node in enumerate(self.model_part.Nodes):
            for k in range(3):
                node_coordinates_setter(node, k, delta)

                self.response_function.Initialize()
                self.response_function.InitializeSolutionStep()

                fd_response_gradient[c * self.block_size + k] = self._CalculateSensitivity(delta)

                node_coordinates_setter(node, k, -delta)

        self._IsVectorRelativelyClose(fd_response_gradient, assembled_adjoint_gradient, 1e-3, 1e-5)

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

    def _AssembleAdjointGradients(self, condition, assembled_vector, element_vector):
            for c, _ in enumerate(condition.GetGeometry()):
                for k in range(3):
                    assembled_vector[self.eq_id_map[condition.Id][c] * self.block_size + k] += element_vector[c * self.block_size + k]

if __name__ == '__main__':
    UnitTest.main()
