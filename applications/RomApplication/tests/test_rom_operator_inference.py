import tempfile

import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.RomApplication.rom_operator_inference import EnergyPreservingNullSpace, OperatorInferenceModel, PhaseInsensitiveScore, QuadraticProducts, SetNodalValues


def DiscreteTrajectory(A, H, B, N, q0, g):
    """q_{k+1} = A q_k + H q2(q_k) + B g_{k+1} + N (g_{k+1} (x) q_k), rows q_0 ... q_K"""
    q = [np.array(q0, dtype=float)]
    for gk in g:
        q.append(A @ q[-1] + H @ QuadraticProducts(q[-1]) + B @ gk + N @ np.kron(gk, q[-1]))
    return np.array(q)


def ImplicitMidpointTrajectory(f, jacobian, q0, g, dt):
    """dq/dt = f(q, g) integrated with the implicit midpoint rule (Newton), rows q_0 ... q_K"""
    q = [np.array(q0, dtype=float)]
    for gk in g:
        x = q[-1].copy()
        for _ in range(50):
            mid = 0.5 * (q[-1] + x)
            dx = np.linalg.solve(np.eye(len(x)) - 0.5 * dt * jacobian(mid, gk), -(x - q[-1] - dt * f(mid, gk)))
            x += dx
            if np.linalg.norm(dx) < 1e-15:
                break
        q.append(x)
    return np.array(q)


def QuadraticJacobian(H, q):
    i, j = np.triu_indices(len(q))
    d = np.zeros((len(i), len(q)))
    d[np.arange(len(i)), i] += q[j]
    d[np.arange(len(i)), j] += q[i]
    return H @ d


class TestRomOperatorInference(KratosUnittest.TestCase):
    """Checks the operator inference ROM (rom_operator_inference.py): exact recovery of known systems, energy preservation, choice of the
    regularization, save and load, and a transient thermal case of the RomApp tests against its FOM."""

    def setUp(self):
        rng = np.random.default_rng(0)
        self.n, self.p = 3, 1
        self.A = 0.8 * np.linalg.qr(rng.standard_normal((self.n, self.n)))[0]
        self.H = 0.05 * rng.standard_normal((self.n, self.n * (self.n + 1) // 2))
        self.B = rng.standard_normal((self.n, self.p))
        self.N = 0.05 * rng.standard_normal((self.n, self.n * self.p))
        steps = 60
        self.inputs = [0.3 * np.sin(0.2 * np.arange(1, steps + 1) * w)[:, None] for w in (1.0, 1.7, 2.3)]
        self.test_input = 0.3 * np.cos(0.25 * np.arange(1, steps + 1))[:, None]

    def testDiscreteRecovery(self):
        # Data of a discrete quadratic system with inputs: the model recovers it (prediction for a new input to round-off)
        trajectories = [DiscreteTrajectory(self.A, self.H, self.B, self.N, np.zeros(self.n), g) for g in self.inputs]
        model = OperatorInferenceModel(("linear", "quadratic", "input", "input_state"), "discrete").Fit(trajectories, self.inputs)
        reference = DiscreteTrajectory(self.A, self.H, self.B, self.N, np.zeros(self.n), self.test_input)
        predicted = model.Predict(reference[0], self.test_input)
        self.assertLess(np.linalg.norm(predicted - reference[1:]) / np.linalg.norm(reference[1:]), 1e-9)

    def testContinuousRecovery(self):
        # Data of dq/dt = A q + H q2 + B g integrated with the implicit midpoint rule: the continuous form recovers it
        dt, A = 0.05, self.A - np.eye(self.n)
        f = lambda q, g: A @ q + self.H @ QuadraticProducts(q) + self.B @ g
        jacobian = lambda q, g: A + QuadraticJacobian(self.H, q)
        trajectories = [ImplicitMidpointTrajectory(f, jacobian, np.zeros(self.n), g, dt) for g in self.inputs]
        model = OperatorInferenceModel(("linear", "quadratic", "input"), "continuous").Fit(trajectories, self.inputs, dt)
        reference = ImplicitMidpointTrajectory(f, jacobian, np.zeros(self.n), self.test_input, dt)
        predicted = model.Predict(reference[0], self.test_input)
        self.assertLess(np.linalg.norm(predicted - reference[1:]) / np.linalg.norm(reference[1:]), 1e-9)

    def testEnergyPreserving(self):
        # Data of dq/dt = H q2 with an energy preserving H: the fitted H satisfies q . H(q (x) q) = 0, and the model conserves |q|
        # with the implicit midpoint rule
        n, dt = 4, 0.05
        rng = np.random.default_rng(1)
        null_space = EnergyPreservingNullSpace(n)
        H = (null_space @ rng.standard_normal(null_space.shape[1])).reshape(n, -1)
        f = lambda q, g: H @ QuadraticProducts(q)
        jacobian = lambda q, g: QuadraticJacobian(H, q)
        trajectories = [ImplicitMidpointTrajectory(f, jacobian, rng.standard_normal(n), np.zeros((80, 0)), dt) for _ in range(3)]
        model = OperatorInferenceModel(("quadratic",), "continuous", energy_preserving=True).Fit(trajectories, time_step=dt)
        for q in rng.standard_normal((10, n)):
            self.assertAlmostEqual(float(q @ (model.operators.T @ QuadraticProducts(q))), 0.0, delta=1e-12 * np.linalg.norm(q) ** 3 * np.abs(model.operators).max())
        q0 = rng.standard_normal(n)
        predicted = model.Predict(q0, number_of_steps=200)
        self.assertLess(np.abs(np.linalg.norm(predicted, axis=1) / np.linalg.norm(q0) - 1.0).max(), 1e-10)

    def testSelectRegularization(self):
        # Noisy data: the candidate with the best leave-one-trajectory-out score is fitted to all the data; if every candidate diverges,
        # the selection fails
        rng = np.random.default_rng(2)
        trajectories = [DiscreteTrajectory(self.A, self.H, self.B, self.N, np.zeros(self.n), g) for g in self.inputs]
        trajectories = [q + 1e-3 * rng.standard_normal(q.shape) for q in trajectories]
        candidates = [(0.0, 0.0), (1e-6, 1e-4), (1e-2, 1e-2)]
        model = OperatorInferenceModel(("linear", "quadratic", "input", "input_state"), "discrete")
        scores = model.SelectRegularization(trajectories, self.inputs, candidates)
        self.assertEqual(set(scores), set(candidates))
        self.assertEqual(scores[model.regularization], min(scores.values()))
        self.assertIsNotNone(model.Predict(trajectories[0][0], self.inputs[0]))
        diverging = OperatorInferenceModel(("linear", "quadratic", "input", "input_state"), "discrete", divergence_bound=1e-6)
        with self.assertRaisesRegex(Exception, "Every candidate regularization diverged"):
            diverging.SelectRegularization(trajectories, self.inputs, candidates)
        self.assertIn("diverged", diverging.failure)

    def testPhaseInsensitiveScore(self):
        # The same oscillation with a phase shift scores (about) zero, a damped one does not
        t = np.arange(400) * 0.05
        reference = np.column_stack([np.sin(2 * np.pi * t), 0.1 * np.cos(2 * np.pi * t)]) + 1.0
        shifted = np.column_stack([np.sin(2 * np.pi * t + 1.0), 0.1 * np.cos(2 * np.pi * t + 1.0)]) + 1.0
        damped = 0.2 * (reference - 1.0) + 1.0
        self.assertLess(PhaseInsensitiveScore(shifted, reference, 0.05), 0.05)
        self.assertGreater(PhaseInsensitiveScore(damped, reference, 0.05), 0.5)

    def testSaveLoad(self):
        trajectories = [DiscreteTrajectory(self.A, self.H, self.B, self.N, np.zeros(self.n), g) for g in self.inputs]
        model = OperatorInferenceModel(("linear", "quadratic", "input", "input_state"), "discrete").Fit(trajectories, self.inputs, regularization=(1e-8, 1e-8))
        with tempfile.TemporaryDirectory() as folder:
            model.Save(folder)
            loaded = OperatorInferenceModel.Load(folder)
        self.assertTrue(np.array_equal(model.Predict(np.zeros(self.n), self.test_input), loaded.Predict(np.zeros(self.n), self.test_input)))

    def testThermalTransient(self):
        # The transient case of thermal_dynamic_test_files (heat flux on top and bottom, convection and radiation on the right) with the
        # heat flux 2000 g(t): FOM snapshots (TEMPERATURE of every node at every step) of four inputs g for the training, and a fifth one
        # (0.75 + 0.5 sin(2 pi t / 1200)) left out. POD of the training snapshots minus the initial state, 6 modes. Both forms predict the
        # left out run within 1.6 % of the FOM (projection on the basis: 0.5 %)
        with KratosUnittest.WorkFolderScope("operator_inference_test_files", __file__):
            data = np.load("thermal_transient_snapshots.npz")
            names = list(data["names"])
            training, test = names[:-1], names[-1]
            time_step = float(data["time_step"])
            u0 = data[f"{training[0]}_states"][0]
            phi = np.linalg.svd(np.vstack([data[f"{name}_states"] - u0 for name in training]).T, full_matrices=False)[0][:, :6]
            trajectories = [(data[f"{name}_states"] - u0) @ phi for name in training]
            inputs = [data[f"{name}_inputs"] for name in training]
            reference = data[f"{test}_states"]
            terms = ("linear", "quadratic", "input", "input_state")
            for form, regularization in (("discrete", (1e-8, 1e-4)), ("continuous", (1e-8, 1e-2))):
                with self.subTest(form=form):
                    model = OperatorInferenceModel(terms, form).Fit(trajectories, inputs, time_step, regularization)
                    q = model.Predict((reference[0] - u0) @ phi, data[f"{test}_inputs"])
                    predicted = u0 + q @ phi.T
                    self.assertLess(np.linalg.norm(predicted - reference[1:]) / np.linalg.norm(reference[1:] - u0), 4e-2)

            # The full order state written to the nodes of the model part
            model_part = KratosMultiphysics.Model().CreateModelPart("Thermal")
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
            KratosMultiphysics.ModelPartIO("../thermal_dynamic_test_files/Square_Radiation_Transient", KratosMultiphysics.IO.READ | KratosMultiphysics.IO.SKIP_TIMER).ReadModelPart(model_part)
            SetNodalValues(model_part, predicted[-1], data["node_ids"], ["TEMPERATURE"])
            for node_id, value in zip(data["node_ids"], predicted[-1]):
                self.assertEqual(model_part.GetNode(int(node_id)).GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), value)


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
