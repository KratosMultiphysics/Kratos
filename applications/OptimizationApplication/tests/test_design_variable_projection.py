import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.opt_projection import CreateProjection

class TestIdentityDesignVariableProjection(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()

        cls.test_model_part =  cls.model.CreateModelPart("test")
        cls.test_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)

        for i in range(10):
            cls.test_model_part.CreateNewNode(i + 1, i, i + 1,0.0)

        cls.position = Kratos.TensorAdaptors.NodePositionTensorAdaptor(cls.test_model_part.Nodes, Kratos.Configuration.Current)
        cls.position.CollectData()

        cls.optimization_problem = OptimizationProblem()
        parameters = Kratos.Parameters("""{
            "type": "identity_projection"
        }""")
        cls.projection = CreateProjection(parameters, cls.optimization_problem)
        cls.projection.SetProjectionSpaces([2.0, 3.0, 4.0], [200.0, 500.0, 1000.0])

    def test_ProjectBackward(self):
        result = self.projection.ProjectBackward(self.position)
        self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 0.0)

    def test_ProjectForward(self):
        result = self.projection.ProjectForward(self.position)
        self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 0.0)

    def test_ForwardProjectionGradient(self):
        result = self.projection.ForwardProjectionGradient(self.position)
        self.assertAlmostEqual(numpy.linalg.norm(result.data), numpy.sqrt(self.test_model_part.NumberOfNodes() * 3))


class TestSigmoidalDesignVariableProjection(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()

        cls.test_model_part =  cls.model.CreateModelPart("test")
        cls.test_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)

        for i in range(10):
            cls.test_model_part.CreateNewNode(i + 1, i * 4, i * 3, i * 2)

        cls.position = Kratos.TensorAdaptors.NodePositionTensorAdaptor(cls.test_model_part.Nodes, Kratos.Configuration.Current)
        cls.position.CollectData()

        cls.optimization_problem = OptimizationProblem()
        parameters = Kratos.Parameters("""{
            "type"          : "sigmoidal_projection",
            "beta_value"    : 2.0,
            "penalty_factor": 1.2
        }""")
        cls.projection = CreateProjection(parameters, cls.optimization_problem)
        cls.projection.SetProjectionSpaces([2.0, 3.0, 4.0], [0.0, 500.0, 2000.0])

    def test_ProjectBackward(self):
        result = self.projection.ProjectBackward(self.position)
        self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 83.57033231673704)

    def test_ProjectForward(self):
        result = self.projection.ProjectForward(self.position)
        self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 9874.416742272824)

    def test_ForwardProjectionGradient(self):
        result = self.projection.ForwardProjectionGradient(self.position)
        self.assertAlmostEqual(numpy.linalg.norm(result.data), 1111.1477428336996)

class TestAdaptiveSigmoidalDesignVariableProjection(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()

        cls.test_model_part =  cls.model.CreateModelPart("test")
        cls.test_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)

        for i in range(10):
            cls.test_model_part.CreateNewNode(i + 1, i * 4, i * 3, i * 2)

        cls.position = Kratos.TensorAdaptors.NodePositionTensorAdaptor(cls.test_model_part.Nodes, Kratos.Configuration.Current)
        cls.position.CollectData()

        cls.parameters = Kratos.Parameters("""{
            "type"          : "adaptive_sigmoidal_projection",
            "initial_value" : 2.0,
            "max_value"     : 3.0,
            "increase_fac"  : 1.2,
            "update_period" : 5,
            "penalty_factor": 1.2
        }""")

    def test_ProjectBackward(self):
        optimization_problem = OptimizationProblem()
        projection = CreateProjection(self.parameters, optimization_problem)
        projection.SetProjectionSpaces([2.0, 3.0, 4.0], [0.0, 500.0, 2000.0])

        for i in range(5):
            optimization_problem.AdvanceStep()
            result = projection.ProjectBackward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 83.57033231673704)

        for i in range(5, 10):
            optimization_problem.AdvanceStep()
            result = projection.ProjectBackward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 82.98662649690897)

        for i in range(10, 15):
            optimization_problem.AdvanceStep()
            result = projection.ProjectBackward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 82.5032449721012)

        for i in range(15, 50):
            optimization_problem.AdvanceStep()
            result = projection.ProjectBackward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 82.40690497800159)

    def test_ProjectForward(self):
        optimization_problem = OptimizationProblem()
        projection = CreateProjection(self.parameters, optimization_problem)
        projection.SetProjectionSpaces([2.0, 3.0, 4.0], [0.0, 500.0, 2000.0])

        for i in range(5):
            optimization_problem.AdvanceStep()
            result = projection.ProjectForward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 9874.416742272824)

        for i in range(5, 10):
            optimization_problem.AdvanceStep()
            result = projection.ProjectForward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 9890.97484074556)

        for i in range(10, 15):
            optimization_problem.AdvanceStep()
            result = projection.ProjectForward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 9905.434754063906)

        for i in range(15, 50):
            optimization_problem.AdvanceStep()
            result = projection.ProjectForward(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 9908.261739693658)

    def test_ForwardProjectionGradient(self):
        optimization_problem = OptimizationProblem()
        projection = CreateProjection(self.parameters, optimization_problem)
        projection.SetProjectionSpaces([2.0, 3.0, 4.0], [0.0, 500.0, 2000.0])

        for i in range(5):
            optimization_problem.AdvanceStep()
            result = projection.ForwardProjectionGradient(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 1108.2457599723648)

        for i in range(5, 10):
            optimization_problem.AdvanceStep()
            result = projection.ForwardProjectionGradient(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 966.0580001392216)

        for i in range(10, 15):
            optimization_problem.AdvanceStep()
            result = projection.ForwardProjectionGradient(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 765.7187019597883)

        for i in range(15, 50):
            optimization_problem.AdvanceStep()
            result = projection.ForwardProjectionGradient(self.position)
            projection.Update()
            self.assertAlmostEqual(numpy.linalg.norm(result.data - self.position.data), 716.3046437643214)

if __name__ == "__main__":
    kratos_unittest.main()