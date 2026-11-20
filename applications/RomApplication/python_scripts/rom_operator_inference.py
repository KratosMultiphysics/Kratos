import itertools
import json
import pathlib

import numpy

import KratosMultiphysics


class OperatorInferenceModel:
    """Reduced model learnt from snapshots (operator inference): a polynomial in the reduced coordinates q and the inputs g, whose
    operators are fitted by regularized least squares. It needs no assembly: the model is integrated online in NumPy, and the full
    order state is u = u0 + Phi q with the basis used to obtain the reduced coordinates of the snapshots.

    form "discrete":   q_{k+1} = f(q_k, g_{k+1}), fitted to consecutive states.
    form "continuous": dq/dt = f(q, g), fitted to (q_{k+1} - q_k) / dt at the midpoint states (q_k + q_{k+1}) / 2 and integrated with the
                       implicit midpoint rule. With energy_preserving the quadratic term satisfies q . H(q (x) q) = 0 (it moves energy
                       between the reduced coordinates and creates none, as the convective term of the Navier-Stokes equations).

    f(q, g) = c + A q + H (q (x) q) + B g + N (g (x) q) + E (g (x) g), with the terms chosen among "constant", "linear", "quadratic",
    "input", "input_state" and "input_quadratic" (q (x) q and g (x) g with their non repeated products only). The state is scaled by one
    factor and each input by its largest absolute value of the training data, so that the regularization weights are comparable between
    problems: lambda_linear for every term but the quadratic one, lambda_quadratic for the quadratic one. For the discrete form without
    energy preservation, the columns of the data matrix are also scaled to unit maximum.
    """

    TERMS = ("constant", "linear", "quadratic", "input", "input_state", "input_quadratic")

    def __init__(self, terms=("linear", "quadratic", "input"), form="discrete", energy_preserving=False, divergence_bound=1.0e6):
        unknown = [term for term in terms if term not in self.TERMS]
        if unknown:
            raise Exception("Unknown terms {}. Available options are: {}.".format(unknown, ", ".join(self.TERMS)))
        if form not in ("discrete", "continuous"):
            raise Exception("Unknown form '{}'. Available options are 'discrete' and 'continuous'.".format(form))
        if energy_preserving and (form != "continuous" or "quadratic" not in terms):
            raise Exception("The energy preserving constraint requires the 'continuous' form and the 'quadratic' term.")
        self.terms = [term for term in self.TERMS if term in terms]
        self.form = form
        self.energy_preserving = energy_preserving
        self.divergence_bound = divergence_bound
        self.operators = None
        self.failure = None

    # ----------------------------------------------------------------------------------------------------------------------------
    # Fit

    def Fit(self, trajectories, inputs=None, time_step=None, regularization=(0.0, 0.0)):
        """Fits the operators to the reduced trajectories (a list of arrays, one per run, rows q_0 ... q_K) and their inputs (a list of
        arrays with one row per step, the input g at the end of each step, K x p). time_step is required by the continuous form.
        regularization: (lambda_linear, lambda_quadratic)."""
        trajectories = [numpy.asarray(q, dtype=numpy.float64) for q in trajectories]
        self.number_of_rom_dofs = trajectories[0].shape[1]
        inputs = self._CheckInputs(trajectories, inputs)
        if self.form == "continuous" and time_step is None:
            raise Exception("The continuous form requires the time_step of the snapshots.")
        self.time_step = time_step
        self.state_scale = max(float(numpy.abs(q).max()) for q in trajectories) or 1.0
        self.input_scale = numpy.maximum(numpy.max([numpy.abs(g).max(axis=0) for g in inputs], axis=0), 1.0e-300) if self.number_of_inputs else numpy.ones(0)

        X, Y = [], []
        for q, g in zip(trajectories, inputs):
            q = q / self.state_scale
            g = g / self.input_scale if self.number_of_inputs else g
            if self.form == "discrete":
                X.append(self._Features(q[:-1], g))
                Y.append(q[1:])
            else:
                X.append(self._Features(0.5 * (q[:-1] + q[1:]), g))
                Y.append((q[1:] - q[:-1]) / time_step)
        X, Y = numpy.vstack(X), numpy.vstack(Y)
        lambda_linear, lambda_quadratic = regularization
        quadratic = self._QuadraticColumns()
        weights = numpy.where(quadratic, lambda_quadratic, lambda_linear) * len(X)

        if self.energy_preserving:
            operators = self._FitEnergyPreserving(X, Y, quadratic, weights)
        else:
            column_scale = numpy.maximum(numpy.abs(X).max(axis=0), 1.0e-300) if self.form == "discrete" else numpy.ones(X.shape[1])
            Xs = X / column_scale
            operators = numpy.linalg.solve(Xs.T @ Xs + numpy.diag(weights), Xs.T @ Y) / column_scale[:, None]
        self.operators = operators
        return self

    def _FitEnergyPreserving(self, X, Y, quadratic, weights):
        """Least squares with q . H(q (x) q) = 0: the non quadratic blocks are eliminated (Schur complement, the same for every output) and
        H is sought in the null space of the constraints on its coefficients."""
        n = self.number_of_rom_dofs
        Xa, Xh = X[:, ~quadratic], X[:, quadratic]
        Ga = Xa.T @ Xa + numpy.diag(weights[~quadratic])
        Gah = Xa.T @ Xh
        Ga_inv_Gah = numpy.linalg.solve(Ga, Gah)
        Ga_inv_XaY = numpy.linalg.solve(Ga, Xa.T @ Y)
        schur = Xh.T @ Xh + numpy.diag(weights[quadratic]) - Gah.T @ Ga_inv_Gah
        rhs = Xh.T @ Y - Gah.T @ Ga_inv_XaY
        null_space = EnergyPreservingNullSpace(n)
        m = Xh.shape[1]
        blocks = null_space.reshape(n, m, -1)
        M = numpy.einsum("lmr,mk,lks->rs", blocks, schur, blocks, optimize=True)
        z = numpy.linalg.solve(M, numpy.einsum("lmr,ml->r", blocks, rhs, optimize=True))
        H = (null_space @ z).reshape(n, m).T
        operators = numpy.zeros((X.shape[1], n))
        operators[quadratic] = H
        operators[~quadratic] = Ga_inv_XaY - Ga_inv_Gah @ H
        return operators

    # ----------------------------------------------------------------------------------------------------------------------------
    # Prediction

    def Predict(self, q0, inputs=None, number_of_steps=None):
        """Reduced trajectory q_1 ... q_K from q0 with the input of each step (K x p, g at the end of each step). Returns None if the
        model diverges (non finite or larger than divergence_bound times the scale of the training states); the reason is in self.failure."""
        if self.operators is None:
            raise Exception("The model must be fitted (Fit) or loaded (Load) before Predict.")
        if self.number_of_inputs:
            inputs = numpy.asarray(inputs, dtype=numpy.float64).reshape(len(inputs), -1) / self.input_scale
            number_of_steps = len(inputs)
        elif number_of_steps is None:
            raise Exception("Without input terms, Predict requires the number_of_steps.")
        q = numpy.asarray(q0, dtype=numpy.float64) / self.state_scale
        self.failure = None
        out = numpy.empty((number_of_steps, self.number_of_rom_dofs))
        for k in range(number_of_steps):
            g = inputs[k:k + 1] if self.number_of_inputs else numpy.zeros((1, 0))
            if self.form == "discrete":
                q = self._Features(q[None], g)[0] @ self.operators
            else:
                q = self._ImplicitMidpointStep(q, g)
                if q is None:
                    self.failure = "the implicit midpoint iteration did not converge at step {}".format(k + 1)
                    return None
            if not numpy.all(numpy.isfinite(q)) or numpy.abs(q).max() > self.divergence_bound:
                self.failure = "the reduced state diverged at step {}".format(k + 1)
                return None
            out[k] = q
        return out * self.state_scale

    def _ImplicitMidpointStep(self, q, g, max_iterations=25, tolerance=1.0e-12):
        dt = self.time_step
        x = q.copy()
        identity = numpy.eye(len(q))
        for _ in range(max_iterations):
            mid = 0.5 * (q + x)
            residual = x - q - dt * (self._Features(mid[None], g)[0] @ self.operators)
            dx = numpy.linalg.solve(identity - 0.5 * dt * self._Jacobian(mid, g[0]), -residual)
            x = x + dx
            if not numpy.all(numpy.isfinite(x)):
                return None
            if numpy.linalg.norm(dx) <= tolerance * max(numpy.linalg.norm(x), 1.0):
                return x
        return None

    def _Jacobian(self, q, g):
        """d f / d q at (q, g), with q and g scaled"""
        n = self.number_of_rom_dofs
        jacobian = numpy.zeros((n, n))
        offset = 0
        for term in self.terms:
            size = self._TermSize(term)
            block = self.operators[offset:offset + size]
            if term == "linear":
                jacobian += block.T
            elif term == "quadratic":
                i, j = numpy.triu_indices(n)
                d = numpy.zeros((len(i), n))
                d[numpy.arange(len(i)), i] += q[j]
                d[numpy.arange(len(i)), j] += q[i]
                jacobian += block.T @ d
            elif term == "input_state":
                jacobian += numpy.einsum("l,lnm->mn", g, block.reshape(self.number_of_inputs, n, n))
            offset += size
        return jacobian

    # ----------------------------------------------------------------------------------------------------------------------------
    # Choice of the regularization

    def SelectRegularization(self, trajectories, inputs=None, candidates=None, time_step=None, score=None):
        """Leave-one-trajectory-out choice of (lambda_linear, lambda_quadratic) among the candidates: for each, the model is fitted to all
        trajectories but one and predicts the one left out from its initial state; the candidate with the lowest mean score is fitted to
        all the trajectories. score(predicted, reference) compares the reduced trajectories q_1 ... q_K (default: relative error); a
        diverged prediction scores infinity. Returns the scores of every candidate."""
        if len(trajectories) < 2:
            raise Exception("SelectRegularization needs at least two trajectories.")
        if candidates is None:
            candidates = list(itertools.product(10.0 ** numpy.arange(-8.0, 1.0, 2.0), 10.0 ** numpy.arange(-8.0, 1.0, 2.0)))
        score = score or RelativeError
        inputs = inputs if inputs is not None else [None] * len(trajectories)
        scores = {}
        for candidate in candidates:
            values = []
            for left in range(len(trajectories)):
                others = [k for k in range(len(trajectories)) if k != left]
                self.Fit([trajectories[k] for k in others], [inputs[k] for k in others] if inputs[0] is not None else None, time_step, candidate)
                q = numpy.asarray(trajectories[left], dtype=numpy.float64)
                predicted = self.Predict(q[0], inputs[left], number_of_steps=len(q) - 1)
                values.append(numpy.inf if predicted is None else float(score(predicted, q[1:])))
            scores[tuple(candidate)] = float(numpy.mean(values))
        best = min(scores, key=scores.get)
        if not numpy.isfinite(scores[best]):
            raise Exception("Every candidate regularization diverged on some left out trajectory.")
        self.Fit(trajectories, inputs if inputs[0] is not None else None, time_step, best)
        self.regularization = best
        return scores

    # ----------------------------------------------------------------------------------------------------------------------------
    # Save and load

    def Save(self, folder, name="OperatorInference"):
        """Saves the model in folder: <name>.npz (operators and scales) and <name>.json (settings)."""
        folder = pathlib.Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        numpy.savez(folder / f"{name}.npz", operators=self.operators, state_scale=self.state_scale, input_scale=self.input_scale)
        settings = {"terms": self.terms, "form": self.form, "energy_preserving": self.energy_preserving, "time_step": self.time_step,
                    "number_of_rom_dofs": self.number_of_rom_dofs, "number_of_inputs": self.number_of_inputs, "divergence_bound": self.divergence_bound}
        with open(folder / f"{name}.json", "w") as settings_file:
            json.dump(settings, settings_file, indent=4)

    @classmethod
    def Load(cls, folder, name="OperatorInference"):
        folder = pathlib.Path(folder)
        with open(folder / f"{name}.json") as settings_file:
            settings = json.load(settings_file)
        model = cls(settings["terms"], settings["form"], settings["energy_preserving"], settings["divergence_bound"])
        data = numpy.load(folder / f"{name}.npz")
        model.operators, model.state_scale, model.input_scale = data["operators"], float(data["state_scale"]), data["input_scale"]
        model.time_step, model.number_of_rom_dofs, model.number_of_inputs = settings["time_step"], settings["number_of_rom_dofs"], settings["number_of_inputs"]
        return model

    # ----------------------------------------------------------------------------------------------------------------------------
    # Features

    def _CheckInputs(self, trajectories, inputs):
        needs_inputs = any(term.startswith("input") for term in self.terms)
        if not needs_inputs:
            self.number_of_inputs = 0
            return [numpy.zeros((len(q) - 1, 0)) for q in trajectories]
        if inputs is None or len(inputs) != len(trajectories):
            raise Exception("The input terms require one input array per trajectory.")
        inputs = [numpy.asarray(g, dtype=numpy.float64).reshape(len(g), -1) for g in inputs]
        for q, g in zip(trajectories, inputs):
            if len(g) != len(q) - 1:
                raise Exception("Each input array must have one row per step ({} rows for a trajectory of {} states).".format(len(q) - 1, len(q)))
        self.number_of_inputs = inputs[0].shape[1]
        return inputs

    def _TermSize(self, term):
        n, p = self.number_of_rom_dofs, self.number_of_inputs
        return {"constant": 1, "linear": n, "quadratic": n * (n + 1) // 2, "input": p, "input_state": p * n, "input_quadratic": p * (p + 1) // 2}[term]

    def _QuadraticColumns(self):
        return numpy.concatenate([numpy.full(self._TermSize(term), term == "quadratic") for term in self.terms])

    def _Features(self, q, g):
        blocks = []
        for term in self.terms:
            if term == "constant":
                blocks.append(numpy.ones((len(q), 1)))
            elif term == "linear":
                blocks.append(q)
            elif term == "quadratic":
                blocks.append(QuadraticProducts(q))
            elif term == "input":
                blocks.append(g)
            elif term == "input_state":
                blocks.append((g[:, :, None] * q[:, None, :]).reshape(len(q), -1))
            elif term == "input_quadratic":
                blocks.append(QuadraticProducts(g))
        return numpy.hstack(blocks)


def QuadraticProducts(x):
    """Products x_i x_j (i <= j) of each row of x"""
    i, j = numpy.triu_indices(x.shape[-1])
    return x[..., i] * x[..., j]


_null_spaces = {}


def EnergyPreservingNullSpace(n):
    """Basis of the coefficients of H (one block of n (n + 1) / 2 per output) such that q . H(q (x) q) = 0 for every q: for each triple
    a <= b <= c the coefficients of q_a q_b q_c in the energy rate cancel."""
    if n not in _null_spaces:
        m = n * (n + 1) // 2
        index = {pair: k for k, pair in enumerate(zip(*numpy.triu_indices(n)))}
        entry = lambda l, i, j: l * m + index[(min(i, j), max(i, j))]
        rows = []
        for a, b, c in itertools.combinations_with_replacement(range(n), 3):
            rows.append(sorted({entry(a, b, c), entry(b, a, c), entry(c, a, b)}))
        C = numpy.zeros((len(rows), n * m))
        for r, columns in enumerate(rows):
            C[r, columns] = 1.0
        _, singular_values, Vt = numpy.linalg.svd(C)
        rank = int((singular_values > 1.0e-10 * singular_values[0]).sum())
        _null_spaces[n] = Vt[rank:].T
    return _null_spaces[n]


def RelativeError(predicted, reference):
    return float(numpy.linalg.norm(predicted - reference) / max(numpy.linalg.norm(reference), 1.0e-300))


def PhaseInsensitiveScore(predicted, reference, time_step, last=None):
    """Score of oscillatory solutions whose phase may drift (e.g. vortex shedding): relative error of the time mean, plus the relative
    errors of the fluctuation (standard deviation) and of the dominant frequency of the component that fluctuates most in the reference,
    over the last 'last' steps (all if None). The pointwise error prefers a damped model to one with the right oscillation and a phase shift."""
    predicted, reference = numpy.asarray(predicted)[-last if last else 0:], numpy.asarray(reference)[-last if last else 0:]
    mean_error = numpy.linalg.norm(predicted.mean(0) - reference.mean(0)) / max(numpy.linalg.norm(reference.mean(0)), 1.0e-300)
    component = int(numpy.argmax(reference.std(0)))
    def Frequency(x):
        x = x - x.mean()
        return numpy.fft.rfftfreq(len(x), time_step)[numpy.argmax(numpy.abs(numpy.fft.rfft(x * numpy.hanning(len(x)))))]
    std_error = abs(predicted[:, component].std() / max(reference[:, component].std(), 1.0e-300) - 1.0)
    frequency_error = abs(Frequency(predicted[:, component]) / max(Frequency(reference[:, component]), 1.0e-300) - 1.0)
    return float(mean_error + std_error + frequency_error)


def SetNodalValues(model_part, values, node_ids, nodal_unknowns, buffer_index=0):
    """Writes a full order vector (e.g. u0 + Phi q) to the nodes of the model part, with the rows ordered as the bases of the RomManager:
    node by node (node_ids, as NodeIds.npy) and, for each node, its nodal_unknowns (as in the RomParameters.json, e.g. ["TEMPERATURE"])."""
    variables = [KratosMultiphysics.KratosGlobals.GetVariable(name) for name in nodal_unknowns]
    values = numpy.asarray(values).reshape(len(node_ids), len(variables))
    for node_id, row in zip(node_ids, values):
        node = model_part.GetNode(int(node_id))
        for variable, value in zip(variables, row):
            node.SetSolutionStepValue(variable, buffer_index, float(value))
