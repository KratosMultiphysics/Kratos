from tkinter import Variable
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from math import isclose

class PerturbationScope:
    def __init__(self, node, variable, delta):
        self.variable = variable
        self.delta = delta
        self.node = node

    def __enter__(self):
        if self.variable == Kratos.SHAPE_SENSITIVITY_X:
            self.node.X += self.delta
        elif self.variable == Kratos.SHAPE_SENSITIVITY_Y:
            self.node.Y += self.delta
        elif self.variable == Kratos.SHAPE_SENSITIVITY_Z:
            self.node.Z += self.delta
        else:
            current_value = self.node.GetSolutionStepValue(self.variable)
            current_value += self.delta
            self.node.SetSolutionStepValue(self.variable, current_value)

    def __exit__(self, exc_type, exc_value, traceback):
        if self.variable == Kratos.SHAPE_SENSITIVITY_X:
            self.node.X -= self.delta
        elif self.variable == Kratos.SHAPE_SENSITIVITY_Y:
            self.node.Y -= self.delta
        elif self.variable == Kratos.SHAPE_SENSITIVITY_Z:
            self.node.Z -= self.delta
        else:
            current_value = self.node.GetSolutionStepValue(self.variable)
            current_value -= self.delta
            self.node.SetSolutionStepValue(self.variable, current_value)

def CalculateSensitivity(perturbed_response_value, ref_value, delta):
    return (perturbed_response_value - ref_value) / delta

def IsVectorRelativelyClose(test_class, vec_a, vec_b, rel_tol, abs_tol):
    test_class.assertEqual(vec_a.Size(), vec_b.Size())
    for i in range(vec_a.Size()):
        if (not isclose(vec_b[i], vec_a[i], rel_tol=rel_tol, abs_tol=abs_tol)):
            msg = "VecA[{:d}] != VecB[{:d}] [ {:f} != {:f} ]. Vectors are: \n\t VecA = {:s}\n\t VecB = {:s}".format(
                i, i, vec_a[i], vec_b[i], str(vec_a), str(vec_b))
            raise AssertionError(msg)

def CheckIsDictRelativelyClose(dict_a, dict_b, rel_tol, abs_tol):
    msg = ""
    if sorted(list(dict_a.keys())) != sorted(list(dict_b.keys())):
        msg = "List of node ids mismatch."
    else:
        for k, v_a in dict_a.items():
            v_b = dict_b[k]
            if len(v_a) != len(v_b):
                msg = "Vector size mismatch at node id {:d} Vectors are: \n\t VecA = {:s}\n\t VecB = {:s}".format(k, str(v_a), str(v_b))
            else:
                for i, v_a_i in enumerate(v_a):
                    if not isclose(v_a_i, v_b[i], rel_tol=rel_tol, abs_tol=abs_tol):
                        msg = "VecA[{:d}] != VecB[{:d}] [ {:f} != {:f} ] at node with id {:d}. Vectors are: \n\t VecA = {:s}\n\t VecB = {:s}".format(
                            i, i, v_a_i, v_b[i], k, str(v_a), str(v_b))

    if msg != "":
        print("Error occured in the comparison of the followings: \n--- Data set A:")
        for node_id in sorted(list(dict_a.keys())):
            print("node {:d} - {:s}".format(node_id, str(dict_a[node_id])))
        print("--- Data set B:")
        for node_id in sorted(list(dict_b.keys())):
            print("node {:d} - {:s}".format(node_id, str(dict_b[node_id])))
        raise AssertionError(msg)

def CalculateAdjointSensitvity(entities_list, element_partial_sensitivty_calculation_lambda):
    sensitivity_dict = {}
    for entity in entities_list:
        number_of_nodes = len(entity.GetGeometry())
        partial_derivatives = element_partial_sensitivty_calculation_lambda(entity)
        number_of_equations = int(partial_derivatives.Size() / number_of_nodes)

        for i, node in enumerate(entity.GetGeometry()):
            if not node.Id in sensitivity_dict.keys():
                sensitivity_dict[node.Id] = [0.0] * number_of_equations

            for equation in range(number_of_equations):
                sensitivity_dict[node.Id][equation] += partial_derivatives[i * number_of_equations + equation]

    return sensitivity_dict

def CalculateFiniteDifferenceSensitivity(nodes_list, variables_list, response_value_calculation_lambda, ref_value, delta):
    sensitivity_dict = {}
    for node in nodes_list:
        sensitivity_dict[node.Id] = [0.0] * len(variables_list)

        for i, sensitivity_variable in enumerate(variables_list):
            if sensitivity_variable is not None:
                with PerturbationScope(node, sensitivity_variable, delta):
                    perturbed_value = response_value_calculation_lambda() / delta
                    sensitivity_dict[node.Id][i] = (perturbed_value - ref_value)

    return sensitivity_dict

def RunResponseSensitivityTest(entities_list, nodes_list, list_of_variables, response_value_calculation_method, element_partial_sensitivty_calculation_method, delta, rel_tol, abs_tol):
    ref_value = response_value_calculation_method() / delta
    adjoint_sensitivities = CalculateAdjointSensitvity(entities_list, element_partial_sensitivty_calculation_method)
    finite_difference_sensitivities = CalculateFiniteDifferenceSensitivity(nodes_list, list_of_variables, response_value_calculation_method, ref_value, delta)
    CheckIsDictRelativelyClose(adjoint_sensitivities, finite_difference_sensitivities, rel_tol, abs_tol)

class TestResidualResponseFunction2D(UnitTest.TestCase):
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
        cls.domain_size = 2
        cls.number_of_nodes = 3
        cls.block_size = 3
        cls.residual_local_size = cls.block_size * cls.number_of_nodes

    @classmethod
    def __CalculateResponseValue(cls):
        return cls.response_function.CalculateValue(cls.model_part)

    def testCalculateFirstDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size)
            self.response_function.CalculateFirstDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Elements,
            self.model_part.Nodes,
            [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.PRESSURE],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-5,
            1e-5,
            1e-5)

    def testCalculateSecondDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size)
            self.response_function.CalculateSecondDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Elements,
            self.model_part.Nodes,
            [Kratos.ACCELERATION_X, Kratos.ACCELERATION_Y, None],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            2e-5,
            1e-7,
            1e-5)

    def testCalculatePartialSensitivity(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.domain_size * self.number_of_nodes, self.domain_size * self.number_of_nodes)
            self.response_function.CalculatePartialSensitivity(
                element,  Kratos.SHAPE_SENSITIVITY, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Elements,
            self.model_part.Nodes,
            [Kratos.SHAPE_SENSITIVITY_X, Kratos.SHAPE_SENSITIVITY_Y],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-7,
            1e-6,
            1e-5)

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
        cls.element = cls.model_part.GetElement(1)
        cls.domain_size = 2
        cls.number_of_nodes = 3
        cls.block_size = 1
        cls.residual_local_size = cls.block_size * cls.number_of_nodes

    @classmethod
    def __CalculateResponseValue(cls):
        return cls.response_function.CalculateValue(cls.model_part)

    def testCalculateFirstDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size)
            self.response_function.CalculateFirstDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Elements,
            self.model_part.Nodes,
            [Kratos.PRESSURE],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-5,
            1e-5,
            1e-5)

    def testCalculateSecondDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size)
            self.response_function.CalculateSecondDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Elements,
            self.model_part.Nodes,
            [None],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            2e-5,
            1e-6,
            1e-5)

    def testCalculatePartialSensitivity(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.domain_size * self.number_of_nodes, self.domain_size * self.number_of_nodes)
            self.response_function.CalculatePartialSensitivity(
                element,  Kratos.SHAPE_SENSITIVITY, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Elements,
            self.model_part.Nodes,
            [Kratos.SHAPE_SENSITIVITY_X, Kratos.SHAPE_SENSITIVITY_Y],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-7,
            1e-6,
            1e-5)


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
        cls.residual_local_size = cls.block_size * cls.number_of_nodes
        cls.nodes_list = {}
        for condition in cls.model_part.Conditions:
            for node in condition.GetGeometry():
                cls.nodes_list[node] = None
        cls.nodes_list = list(cls.nodes_list.keys())

    @classmethod
    def __CalculateResponseValue(cls):
        return cls.response_function.CalculateValue(cls.model_part)

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

        self.assertAlmostEqual(value / total_area, self.__CalculateResponseValue(), 9)

    def testCalculateFirstDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size)
            self.response_function.CalculateFirstDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Conditions,
            self.nodes_list,
            [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, None],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-8,
            1e-5,
            1e-5)

    def testCalculateSecondDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size)
            self.response_function.CalculateSecondDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Conditions,
            self.nodes_list,
            [None, None, None],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            2e-5,
            1e-6,
            1e-5)

    def testCalculatePartialSensitivity(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.domain_size * self.number_of_nodes, self.domain_size * self.number_of_nodes)
            self.response_function.CalculatePartialSensitivity(
                element,  Kratos.SHAPE_SENSITIVITY, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.model_part.Conditions,
            self.nodes_list,
            [Kratos.SHAPE_SENSITIVITY_X, Kratos.SHAPE_SENSITIVITY_Y, Kratos.SHAPE_SENSITIVITY_Z],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-5,
            1e-3,
            1e-5)

class TestDragResponseFunction2DSteady(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.domain_size = 2
        cls.number_of_nodes = 3
        cls.block_size = 3
        cls.residual_local_size = cls.block_size * cls.number_of_nodes

        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.ProcessInfo[Kratos.TIME] = 0.1
        cls.model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.04
        cls.model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 0.3
        cls.model_part.ProcessInfo[Kratos.OSS_SWITCH] = 0
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = cls.domain_size

        cls.model_part.SetBufferSize(1)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ADVPROJ)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_1)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_2)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_3)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_SCALAR_1)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,cls.model_part)

        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_X, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_SCALAR_1, cls.model_part)

        for node in cls.model_part.Nodes:
            node.Fix(Kratos.VELOCITY_X)
            node.Fix(Kratos.VELOCITY_Y)
            node.Fix(Kratos.VELOCITY_Z)
            node.Fix(Kratos.PRESSURE)

        sub_model_part = cls.model_part.CreateSubModelPart("response_surface")
        sub_model_part.AddNodes([1, 2, 3, 4])

        prop = cls.model_part.GetProperties()[0]
        prop[Kratos.DENSITY] = 1.5
        prop[Kratos.DYNAMIC_VISCOSITY] = 1.2e-4
        prop[Kratos.CONSTITUTIVE_LAW] = KratosCFD.Newtonian2DLaw()

        cls.primal_model_part = cls.model.CreateModelPart("test_primal")
        cls.primal_model_part.ProcessInfo[Kratos.TIME] = 0.1
        cls.primal_model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.04
        cls.primal_model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 0.3
        cls.primal_model_part.ProcessInfo[Kratos.OSS_SWITCH] = 0
        cls.primal_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = cls.domain_size
        for node in cls.model_part.Nodes:
            cls.primal_model_part.AddNode(node, 0)
        cls.primal_model_part.CreateNewElement("QSVMS2D3N", 1, [2, 3, 1], prop)
        cls.primal_model_part.CreateNewElement("QSVMS2D3N", 2, [3, 4, 1], prop)

        cls.adjoint_model_part = cls.model.CreateModelPart("test_adjoint")
        cls.adjoint_model_part.ProcessInfo[Kratos.TIME] = 0.1
        cls.adjoint_model_part.ProcessInfo[Kratos.DELTA_TIME] = -0.04
        cls.adjoint_model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 0.3
        cls.adjoint_model_part.ProcessInfo[Kratos.OSS_SWITCH] = 0
        cls.adjoint_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = cls.domain_size
        for node in cls.model_part.Nodes:
            cls.adjoint_model_part.AddNode(node, 0)
        cls.adjoint_model_part.CreateNewElement("QSVMSAdjoint2D3N", 3, [2, 3, 1], prop)
        cls.adjoint_model_part.CreateNewElement("QSVMSAdjoint2D3N", 4, [3, 4, 1], prop)

        sub_model_part = cls.adjoint_model_part.CreateSubModelPart("response_surface")
        sub_model_part.AddNodes([1, 2, 3, 4])

        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VELOCITY, 1.0, 10.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.PRESSURE, 0.0, 10.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.MESH_VELOCITY, 0.0, 1.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.ACCELERATION, 0.0, 0.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.DENSITY, 1.5, 1.5, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VISCOSITY, 1.2e-4, 1.2e-4, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.BODY_FORCE, 10.0, 10.0, 0)

        cls.primal_scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(1e-9, 1e-9, 2)
        cls.primal_scheme.Initialize(cls.primal_model_part)
        cls.primal_scheme.InitializeElements(cls.primal_model_part)
        cls.primal_scheme.InitializeSolutionStep(cls.primal_model_part, Kratos.CompressedMatrix(), Kratos.Vector(), Kratos.Vector())

        cls.response_function = KratosCFD.DragResponseFunction2D(Kratos.Parameters("""{
            "structure_model_part_name": "response_surface",
            "drag_direction"           : [1.0, 0.0, 0.0],
            "start_time"               : 0.0
        }"""),  cls.model_part)
        cls.response_function.Initialize()

        cls.adjoint_scheme = KratosCFD.SimpleSteadyAdjointScheme(cls.response_function, cls.domain_size, cls.block_size)
        cls.adjoint_scheme.Initialize(cls.adjoint_model_part)
        cls.adjoint_scheme.InitializeElements(cls.adjoint_model_part)
        cls.adjoint_scheme.InitializeSolutionStep(cls.adjoint_model_part, Kratos.CompressedMatrix(), Kratos.Vector(), Kratos.Vector())

    @classmethod
    def __CalculateResponseValue(cls):
        cls.primal_scheme.FinalizeSolutionStep(cls.primal_model_part, Kratos.CompressedMatrix(), Kratos.Vector(), Kratos.Vector())
        return cls.response_function.CalculateValue(cls.primal_model_part)

    def testCalculateFirstDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix()
            element.CalculateFirstDerivativesLHS(residual_gradient, self.adjoint_model_part.ProcessInfo)
            self.response_function.CalculateFirstDerivativesGradient(
                element, residual_gradient * -1.0, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.adjoint_model_part.Elements,
            self.primal_model_part.Nodes,
            [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.PRESSURE],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-6,
            1e-10,
            1e-6)

    def testCalculateSecondDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size, 0.0)
            self.response_function.CalculateSecondDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.adjoint_model_part.Elements,
            self.primal_model_part.Nodes,
            [Kratos.ACCELERATION_X, Kratos.ACCELERATION_Y, None],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-5,
            1e-6,
            1e-5)

    def testCalculatePartialSensitivity(self):
        delta = 1e-6

        ref_value = self.__CalculateResponseValue() / delta
        finite_difference_sensitivities = CalculateFiniteDifferenceSensitivity(self.primal_model_part.Nodes, [Kratos.SHAPE_SENSITIVITY_X, Kratos.SHAPE_SENSITIVITY_Y], lambda : self.__CalculateResponseValue(), ref_value, delta)

        sensitivity_builder_settings = Kratos.Parameters("""{
            "sensitivity_model_part_name": "response_surface",
            "nodal_solution_step_sensitivity_variables": [
                "SHAPE_SENSITIVITY"
            ],
            "build_mode": "sum",
            "nodal_solution_step_sensitivity_calculation_is_thread_safe": true
        }""")
        sensitivity_builder_scheme = KratosCFD.SimpleSteadySensitivityBuilderScheme(self.domain_size, self.block_size)
        sensitivity_builder = Kratos.SensitivityBuilder(
            sensitivity_builder_settings,
            self.adjoint_model_part,
            self.response_function,
            sensitivity_builder_scheme)
        sensitivity_builder.Initialize()
        sensitivity_builder.InitializeSolutionStep()
        sensitivity_builder.UpdateSensitivities()

        adjoint_sensitivities = {}
        for node in self.adjoint_model_part.Nodes:
            sensitivities = node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)
            adjoint_sensitivities[node.Id] = [-sensitivities[0], -sensitivities[1]]

        CheckIsDictRelativelyClose(adjoint_sensitivities, finite_difference_sensitivities, 1e-5, 1e-6)

class TestDragResponseFunction2DBossak(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.domain_size = 2
        cls.number_of_nodes = 3
        cls.block_size = 3
        cls.residual_local_size = cls.block_size * cls.number_of_nodes

        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.model_part.ProcessInfo[Kratos.TIME] = 0.1
        cls.model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.04
        cls.model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 0.3
        cls.model_part.ProcessInfo[Kratos.OSS_SWITCH] = 0
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = cls.domain_size

        cls.model_part.SetBufferSize(2)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ADVPROJ)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_1)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_2)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_3)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_SCALAR_1)
        cls.model_part.AddNodalSolutionStepVariable(KratosCFD.AUX_ADJOINT_FLUID_VECTOR_1)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,cls.model_part)

        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_X, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_SCALAR_1, cls.model_part)

        for node in cls.model_part.Nodes:
            node.Fix(Kratos.VELOCITY_X)
            node.Fix(Kratos.VELOCITY_Y)
            node.Fix(Kratos.VELOCITY_Z)
            node.Fix(Kratos.PRESSURE)

        sub_model_part = cls.model_part.CreateSubModelPart("response_surface")
        sub_model_part.AddNodes([1, 2, 3, 4])

        prop = cls.model_part.GetProperties()[0]
        prop[Kratos.DENSITY] = 1.5
        prop[Kratos.DYNAMIC_VISCOSITY] = 1.2e-4
        prop[Kratos.CONSTITUTIVE_LAW] = KratosCFD.Newtonian2DLaw()

        cls.primal_model_part = cls.model.CreateModelPart("test_primal")
        cls.primal_model_part.ProcessInfo[Kratos.TIME] = 0.1
        cls.primal_model_part.ProcessInfo[Kratos.DELTA_TIME] = 0.04
        cls.primal_model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 0.3
        cls.primal_model_part.ProcessInfo[Kratos.OSS_SWITCH] = 0
        cls.primal_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = cls.domain_size
        for node in cls.model_part.Nodes:
            cls.primal_model_part.AddNode(node, 0)
        cls.primal_model_part.CreateNewElement("QSVMS2D3N", 1, [2, 3, 1], prop)
        cls.primal_model_part.CreateNewElement("QSVMS2D3N", 2, [3, 4, 1], prop)

        cls.adjoint_model_part = cls.model.CreateModelPart("test_adjoint")
        cls.adjoint_model_part.ProcessInfo[Kratos.TIME] = 0.1
        cls.adjoint_model_part.ProcessInfo[Kratos.DELTA_TIME] = -0.04
        cls.adjoint_model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = 0.3
        cls.adjoint_model_part.ProcessInfo[Kratos.OSS_SWITCH] = 0
        cls.adjoint_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = cls.domain_size
        for node in cls.model_part.Nodes:
            cls.adjoint_model_part.AddNode(node, 0)
        cls.adjoint_model_part.CreateNewElement("QSVMSAdjoint2D3N", 3, [2, 3, 1], prop)
        cls.adjoint_model_part.CreateNewElement("QSVMSAdjoint2D3N", 4, [3, 4, 1], prop)

        sub_model_part = cls.adjoint_model_part.CreateSubModelPart("response_surface")
        sub_model_part.AddNodes([1, 2, 3, 4])

        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VELOCITY, 1.0, 10.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.PRESSURE, 0.0, 10.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.MESH_VELOCITY, 0.0, 1.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.ACCELERATION, 0.0, 1.0, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.DENSITY, 1.5, 1.5, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VISCOSITY, 1.2e-4, 1.2e-4, 0)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.BODY_FORCE, 10.0, 10.0, 0)

        cls.adjoint_model_part.CloneTimeStep(2.0)
        cls.adjoint_model_part.CloneTimeStep(1.96)

        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VELOCITY, 1.0, 10.0, 1)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.PRESSURE, 0.0, 10.0, 1)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.MESH_VELOCITY, 0.0, 1.0, 1)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.ACCELERATION, 0.0, 1.0, 1)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.DENSITY, 1.5, 1.5, 1)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.VISCOSITY, 1.2e-4, 1.2e-4, 1)
        KratosCFD.FluidTestUtilities.RandomFillNodalHistoricalVariable(cls.model_part, Kratos.BODY_FORCE, 10.0, 10.0, 1)

        cls.primal_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(-0.3, 0, cls.domain_size)
        cls.primal_scheme.Initialize(cls.primal_model_part)
        cls.primal_scheme.InitializeElements(cls.primal_model_part)
        cls.primal_scheme.InitializeSolutionStep(cls.primal_model_part, Kratos.CompressedMatrix(), Kratos.Vector(), Kratos.Vector())

        cls.response_function = KratosCFD.DragResponseFunction2D(Kratos.Parameters("""{
            "structure_model_part_name": "response_surface",
            "drag_direction"           : [1.0, 0.0, 0.0],
            "start_time"               : 0.0
        }"""),  cls.model_part)
        cls.response_function.Initialize()

        cls.adjoint_scheme = KratosCFD.VelocityBossakAdjointScheme(Kratos.Parameters("""{"alpha_bossak": -0.3}"""), cls.response_function, cls.domain_size, cls.block_size)
        cls.adjoint_scheme.Initialize(cls.adjoint_model_part)
        cls.adjoint_scheme.InitializeElements(cls.adjoint_model_part)
        cls.adjoint_scheme.InitializeSolutionStep(cls.adjoint_model_part, Kratos.CompressedMatrix(), Kratos.Vector(), Kratos.Vector())

    @classmethod
    def __CalculateResponseValue(cls):
        cls.primal_scheme.FinalizeSolutionStep(cls.primal_model_part, Kratos.CompressedMatrix(), Kratos.Vector(), Kratos.Vector())
        return cls.response_function.CalculateValue(cls.primal_model_part)

    def testCalculateFirstDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix()
            element.CalculateFirstDerivativesLHS(residual_gradient, self.adjoint_model_part.ProcessInfo)
            self.response_function.CalculateFirstDerivativesGradient(
                element, residual_gradient * -1.0, response_gradient, self.model_part.ProcessInfo)
            return response_gradient

        RunResponseSensitivityTest(
            self.adjoint_model_part.Elements,
            self.primal_model_part.Nodes,
            [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.PRESSURE],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-6,
            1e-5,
            1e-6)

    def testCalculateSecondDerivativesGradient(self):
        def element_partial_sensitivty_calculation_method(element):
            response_gradient = Kratos.Vector()
            residual_gradient = Kratos.Matrix(self.residual_local_size, self.residual_local_size)
            element.CalculateSecondDerivativesLHS(residual_gradient, self.adjoint_model_part.ProcessInfo)
            self.response_function.CalculateSecondDerivativesGradient(
                element, residual_gradient, response_gradient, self.model_part.ProcessInfo)
            return response_gradient * -1.3

        RunResponseSensitivityTest(
            self.adjoint_model_part.Elements,
            self.primal_model_part.Nodes,
            [Kratos.ACCELERATION_X, Kratos.ACCELERATION_Y, None],
            lambda : self.__CalculateResponseValue(),
            element_partial_sensitivty_calculation_method,
            1e-5,
            1e-6,
            1e-5)

    def testCalculatePartialSensitivity(self):
        delta = 1e-6

        ref_value = self.__CalculateResponseValue() / delta
        finite_difference_sensitivities = CalculateFiniteDifferenceSensitivity(self.primal_model_part.Nodes, [Kratos.SHAPE_SENSITIVITY_X, Kratos.SHAPE_SENSITIVITY_Y], lambda : self.__CalculateResponseValue(), ref_value, delta)

        sensitivity_builder_settings = Kratos.Parameters("""{
            "sensitivity_model_part_name": "response_surface",
            "nodal_solution_step_sensitivity_variables": [
                "SHAPE_SENSITIVITY"
            ],
            "build_mode": "sum",
            "nodal_solution_step_sensitivity_calculation_is_thread_safe": true
        }""")
        sensitivity_builder_scheme = KratosCFD.VelocityBossakSensitivityBuilderScheme(-0.3, self.domain_size, self.block_size)
        sensitivity_builder = Kratos.SensitivityBuilder(
            sensitivity_builder_settings,
            self.adjoint_model_part,
            self.response_function,
            sensitivity_builder_scheme)
        sensitivity_builder.Initialize()
        sensitivity_builder.InitializeSolutionStep()
        sensitivity_builder.UpdateSensitivities()

        adjoint_sensitivities = {}
        for node in self.adjoint_model_part.Nodes:
            sensitivities = node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)
            adjoint_sensitivities[node.Id] = [-sensitivities[0], -sensitivities[1]]

        CheckIsDictRelativelyClose(adjoint_sensitivities, finite_difference_sensitivities, 1e-5, 1e-6)

if __name__ == '__main__':
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    UnitTest.main()
