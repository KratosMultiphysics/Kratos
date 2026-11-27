# PROTOTYPE: the Future ROM as the solution strategy of a legacy Kratos solver, with no builder and solver.
#
#   - The solver (any PythonSolver with _CreateSolutionStrategy / _GetScheme / _GetConvergenceCriterion) is the FOM one, untouched: its
#     own step methods run and call this strategy. Coupled solvers get one strategy per sub-solver (UseFutureRomStrategies).
#   - Future stack: DOFs, system and assembly by a Future scheme (LegacySchemeAdapter, C++), projection and reduced solve by the C++
#     RomProjector of the RomApp, decoder of the ROM folder (CreateDecoder: any registered decoder) on the DOFs of this solver.
#   - Physics taken from the solver AS OBJECTS: its time scheme, wrapped by the adapter (contribution of each entity, Predict, Update,
#     operations of each step and iteration), and its convergence criterion.
#   - The Newton loop is here, in Python, in the sequence of the legacy ResidualBasedNewtonRaphsonStrategy.
import contextlib
import io
import json
import os
import pathlib

import numpy
import scipy.sparse

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics import scipy_conversion_tools
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.RomApplication.future.rom_decoders import CreateDecoder, GetBasisRows
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition


# --- Future scheme of a solver: native when there is one, else the adapter of its legacy scheme ---------------------------------------
# Legacy scheme type -> (class of KratosROM.Future for the ROM runs, class for the HROM runs) built with (model_part, settings).
# A scheme not registered here is wrapped by the LegacySchemeAdapter, so every solver works now; registering the Future twin of a
# scheme when the core provides it is the whole switch.
_NATIVE_SCHEMES = {}


def RegisterNativeScheme(legacy_scheme_name, rom_scheme_class_name, hrom_scheme_class_name):
    _NATIVE_SCHEMES[legacy_scheme_name] = (rom_scheme_class_name, hrom_scheme_class_name)


# The legacy static scheme and the Future static scheme do the same: entities alone, u += dx
RegisterNativeScheme("ResidualBasedIncrementalUpdateStaticScheme", "RomStaticScheme", "HRomStaticScheme")


def UnregisterNativeScheme(legacy_scheme_name):
    _NATIVE_SCHEMES.pop(legacy_scheme_name, None)


def NativeSchemes():
    """The legacy scheme types that have a native Future scheme registered: {legacy name: (ROM class, HROM class)}"""
    return dict(_NATIVE_SCHEMES)


def CreateFutureScheme(model_part, settings, legacy_scheme, run_hrom=False, scheme="auto"):
    """Returns the Future scheme and how it was obtained ("native", "adapter" or "given").
    scheme: "auto" (the native scheme registered for the legacy one, else the adapter), "native" (error if none is registered),
    "adapter" (always the adapter of the legacy scheme), or a function (model_part, settings, legacy_scheme, run_hrom) -> Future scheme."""
    legacy_scheme_name = type(legacy_scheme).__name__
    if callable(scheme):
        return scheme(model_part, settings, legacy_scheme, run_hrom), "given"
    if scheme not in ("auto", "native", "adapter"):
        raise Exception("Unknown 'scheme' \'{}\'. Options: 'auto', 'native', 'adapter' or a function.".format(scheme))
    if scheme == "native" and legacy_scheme_name not in _NATIVE_SCHEMES:
        raise Exception("No native Future scheme is registered for \'{}\'. Registered: {}. Use 'auto' or 'adapter', or register one "
                        "with RegisterNativeScheme.".format(legacy_scheme_name, ", ".join(_NATIVE_SCHEMES) or "none"))
    if scheme != "adapter" and legacy_scheme_name in _NATIVE_SCHEMES:
        return getattr(KratosROM.Future, _NATIVE_SCHEMES[legacy_scheme_name][int(run_hrom)])(model_part, settings), "native"
    adapter_class = KratosROM.Future.HRomLegacySchemeAdapter if run_hrom else KratosROM.Future.RomLegacySchemeAdapter
    return adapter_class(model_part, settings, legacy_scheme), "adapter"


class FutureRomStrategy:
    def __init__(self, solver, rom_data_folder="rom_data", criterion_residual="zero", run_hrom=False, hrom_weights_set=0, train_hrom=False,
                 use_native_schemes=True, convergence_criterion="solver", residual_relative_tolerance=1.0e-9, residual_absolute_tolerance=1.0e-9,
                 update="increment", scheme="auto", decoder_factory=None, name=None, snapshots=None, initial_reduced_coordinates=None,
                 projection="auto", lspg_solving_technique=None, monotonicity_preserving=None):
        """criterion_residual: the residual vector given to the convergence criterion. "zero" as the legacy ROM (its builder never
        builds the full residual, so a residual criterion converges at once); "full" the assembled FOM residual.
        run_hrom: only the HROM entities of the ROM folder are assembled, scaled by the column hrom_weights_set of their weights (one
        column per physics in coupled problems).
        train_hrom: the residual of each entity, through the scheme, is projected on the Jacobian of the decoder at the end of each
        step (hrom_residual_snapshots: the residual snapshots of the HROM training).
        scheme: which Future scheme stands for the scheme of the solver: "auto", "native", "adapter" or a function (see
        CreateFutureScheme). use_native_schemes=False is the same as scheme="adapter".
        decoder_factory: function (rom_data_folder, dof_set) -> decoder, instead of CreateDecoder (the decoder of the ROM folder).
        name: name of the solver in a coupled one (e.g. "fluid_solver"), for the reports.
        snapshots: full order states to follow instead of solving (e.g. of the FOM): matrix with one column per step and the rows of
        the bases of the ROM folder (nodes of NodeIds.npy times 'nodal_unknowns'). In each step the state of this solver is set to
        the snapshot projected on its decoder, D(Encode(u)), through the update of the scheme, with no reduced system. With train_hrom
        this gives the HROM training residuals without any ROM run.
        convergence_criterion: "solver", the criterion object of the solver (as the legacy ROM), or "reduced_residual_criterion", the
        one of the Future ROM (#14955): norm of the projected residual before each solve, relative to the one of the first iteration
        of the step (residual_relative_tolerance) or absolute, divided by the square root of the number of ROM DOFs.
        update: how the new reduced coordinates reach the DOFs. "increment": D(q + dq) - D(q) added by the scheme on top of the current
        state (the legacy ROM with a linear decoder; keeps the predictor of the scheme and any part of the state outside the decoder).
        "absolute": the state is set to D(q) before each build (the FutureRomSolver of #14850; the predictor is overwritten).
        projection: the test basis W of the reduced system W^T A Phi dq = W^T b, with Phi the Jacobian of the decoder. "galerkin" (W = Phi),
        "petrov_galerkin" (W = Psi, the LeftBasisMatrix.npy of the ROM folder with its 'petrov_galerkin_number_of_rom_dofs' columns; solved
        by least squares when it has more columns than Phi) or "lspg" (W = A Phi: least squares of A Phi dq = b). "auto": the
        'projection_strategy' of the RomParameters.json. lspg_solving_technique: "qr_decomposition" (least squares on A Phi) or
        "normal_equations" ((A Phi)^T A Phi); None: the one of the 'rom_bns_settings' of the RomParameters.json, else "normal_equations".
        monotonicity_preserving: the correction of the legacy ROM builders on the assembled system before it is projected: every positive
        off-diagonal entry is removed and added to the diagonal (and the residual changed accordingly), which makes the matrix an
        M-matrix. None: the 'monotonicity_preserving' of the 'rom_bns_settings' of the RomParameters.json, else False.
        initial_reduced_coordinates: None (q = 0 and the initial state as it is, enough for a linear decoder, whose Jacobian does not
        depend on q), "encode" (the q whose decoding is closest to the initial state on the free DOFs, by Gauss-Newton) or a function
        (decoder, dof_set) -> q. With the last two the free DOFs are set to D(q) before the first step, as the legacy ANN and RBF
        enhanced ROMs do: a nonlinear decoder needs the reduced coordinates of the state it starts from."""
        self.solver = solver
        self.model_part = solver.GetComputingModelPart()
        self.rom_data_folder = rom_data_folder
        self.criterion_residual = criterion_residual
        self.criterion = solver._GetConvergenceCriterion()
        settings = solver.settings
        self.max_iterations = next(settings[name].GetInt() for name in ("maximum_iterations", "max_iteration") if settings.Has(name))

        # Future scheme with the time integration of the scheme of the solver
        scheme_settings = KratosMultiphysics.Parameters("""{
            "name" : "static_scheme",
            "build_settings" : {"name" : "block_builder"},
            "echo_level" : 0,
            "move_mesh" : false,
            "reform_dofs_at_each_step" : false
        }""")
        if settings.Has("move_mesh_flag"):
            scheme_settings["move_mesh"].SetBool(settings["move_mesh_flag"].GetBool())
        self.legacy_scheme = solver._GetScheme()
        self.run_hrom, self.hrom_weights_set = run_hrom, hrom_weights_set
        if train_hrom and run_hrom:
            raise Exception("The HROM training requires the ROM, not the HROM.")
        if convergence_criterion not in ("solver", "reduced_residual_criterion") or update not in ("increment", "absolute"):
            raise Exception("Unknown 'convergence_criterion' or 'update'.")
        self.convergence_criterion, self.update = convergence_criterion, update
        self.residual_relative_tolerance, self.residual_absolute_tolerance = residual_relative_tolerance, residual_absolute_tolerance
        self.train_hrom, self.hrom_residual_snapshots = train_hrom, []
        self.name, self.decoder_factory = name, decoder_factory
        self.snapshots = snapshots
        with open(pathlib.Path(rom_data_folder) / "RomParameters.json") as parameter_file:
            rom_parameters = json.load(parameter_file)
        if projection == "auto":
            projection = {"lspg": "lspg", "petrov_galerkin": "petrov_galerkin"}.get(rom_parameters.get("projection_strategy", "galerkin"), "galerkin")
        if projection not in ("galerkin", "petrov_galerkin", "lspg"):
            raise Exception("Unknown 'projection' \'{}\'. Options: 'auto', 'galerkin', 'petrov_galerkin', 'lspg'.".format(projection))
        self.projection = projection
        self.lspg_solving_technique = lspg_solving_technique or rom_parameters["rom_settings"].get("rom_bns_settings", {}).get("solving_technique", "normal_equations")
        if self.lspg_solving_technique not in ("qr_decomposition", "normal_equations"):
            raise Exception("Unknown 'lspg_solving_technique' \'{}\'.".format(self.lspg_solving_technique))
        if projection == "lspg" and train_hrom:
            raise Exception("The HROM training is not available with the LSPG projection.")
        if projection == "lspg" and run_hrom and self.lspg_solving_technique != "normal_equations":
            raise Exception("The LSPG HROM requires the 'normal_equations' solving technique (as the legacy one).")
        self.monotonicity_preserving = bool(rom_parameters["rom_settings"].get("rom_bns_settings", {}).get("monotonicity_preserving", False)) if monotonicity_preserving is None else monotonicity_preserving
        self.petrov_galerkin_number_of_rom_dofs = rom_parameters["rom_settings"].get("petrov_galerkin_number_of_rom_dofs", 0)
        self.initial_reduced_coordinates, self.initial_state_encoded = initial_reduced_coordinates, False
        self.scheme, self.scheme_origin = CreateFutureScheme(self.model_part, scheme_settings, self.legacy_scheme, run_hrom,
                                                             scheme if use_native_schemes else "adapter")
        self.strategy_data = KratosMultiphysics.Future.ImplicitStrategyData()
        self.projector = KratosROM.Future.RomProjector(self.scheme, self.strategy_data)

        # Arguments of the legacy interface of the criterion (the matrix is not used)
        self.space = KratosMultiphysics.SparseSpace()
        self.criterion_A, self.criterion_dx, self.criterion_b = self.space.CreateEmptyMatrixPointer(), self.space.CreateEmptyVectorPointer(), self.space.CreateEmptyVectorPointer()
        self.initialized = False
        self.decoder = None
        self.iterations = []  # Newton iterations of each step
        self.local_model_steps = []  # model of each step, for a decoder made of local models
        self.rom_solutions = []

    # --- interface the solvers call on their strategy -------------------------------------------------------------------------------
    def Initialize(self):
        if not self.initialized:
            self.scheme.Initialize(self.strategy_data)  # DOFs, equation ids, system; scheme, elements and conditions
            self.criterion.Initialize(self.model_part)
            self.dof_set = self.strategy_data.GetDofSet()
            for vector in (self.criterion_dx, self.criterion_b):
                self.space.ResizeVector(vector, len(self.dof_set))
                self.space.SetToZeroVector(vector)
            # Decoder of the ROM folder with its rows on the (effective) DOFs of this solver
            self.decoder = (self.decoder_factory or CreateDecoder)(self.rom_data_folder, self.strategy_data.GetEffectiveDofSet())
            self.q = numpy.zeros(self.decoder.NumberOfRomDofs())
            if self.projection == "petrov_galerkin":
                left_modes = numpy.load(pathlib.Path(self.rom_data_folder) / "LeftBasisMatrix.npy")
                rows = GetBasisRows(self.rom_data_folder, self.strategy_data.GetEffectiveDofSet())
                self.left_basis = numpy.ascontiguousarray(left_modes[rows, :self.petrov_galerkin_number_of_rom_dofs or left_modes.shape[1]])
            # Distance of the decoder at q = 0 to the initial state (with "increment" it stays in the state as a constant offset)
            self.initial_decoder_offset = float(numpy.linalg.norm(self.decoder.Decode(self.q) - numpy.array(self.strategy_data.GetEffectiveDofSet().GetValues())))
            if self.run_hrom:
                self._SetHRomEntities()
            if self.snapshots is not None:
                self.snapshot_rows = GetBasisRows(self.rom_data_folder, self.strategy_data.GetEffectiveDofSet())
            self.initialized = True

    def _EncodeInitialState(self):
        """Reduced coordinates of the initial state and free DOFs on the decoder, once, before the first step"""
        if self.initial_state_encoded or self.initial_reduced_coordinates is None:
            return
        self.initial_state_encoded = True
        dof_set = self.strategy_data.GetEffectiveDofSet()
        state = numpy.array(dof_set.GetValues())
        free = numpy.array([not dof.IsFixed() for dof in dof_set])
        if callable(self.initial_reduced_coordinates):
            q = numpy.array(self.initial_reduced_coordinates(self.decoder, dof_set), dtype=numpy.float64)
        elif self.initial_reduced_coordinates == "encode":
            q = numpy.zeros(self.decoder.NumberOfRomDofs())
            for _ in range(20):
                dq = numpy.linalg.lstsq(self.decoder.Jacobian(q)[free], (state - self.decoder.Decode(q))[free], rcond=None)[0]
                q = q + dq
                if numpy.linalg.norm(dq) <= 1.0e-12 * max(numpy.linalg.norm(q), 1.0e-300):
                    break
        else:
            raise Exception("Unknown 'initial_reduced_coordinates'. Options: None, 'encode' or a function (decoder, dof_set) -> q.")
        self.q = q
        state[free] = self.decoder.Decode(q)[free]
        dof_set.SetValues(state)

    def InitializeSolutionStep(self):
        self.Initialize()
        self._EncodeInitialState()
        self.scheme.InitializeSolutionStep(self.strategy_data)
        self.criterion.InitializeSolutionStep(self.model_part, self.dof_set, self.criterion_A, self.criterion_dx, self._CriterionResidual())
        self.decoder.InitializeSolutionStep()
        # A single decoder keeps the reduced coordinates; one made of local models (LocalRomDecoder) may switch model here and encode
        # the current state in the new one
        if hasattr(self.decoder, "UpdateReducedCoordinates"):
            model = self.decoder.GetActiveModel() if self.local_model_steps else None
            self.q = numpy.array(self.decoder.UpdateReducedCoordinates(self.q, self.strategy_data.GetEffectiveDofSet(), self.model_part.ProcessInfo), dtype=numpy.float64)
            if hasattr(self.decoder, "GetActiveModel"):
                # With "increment" the state is put on the new model at a switch (and at the first step): the increments that follow
                # are the ones of that model, and what the old model left out of it would stay in the state otherwise
                if self.update == "increment" and self.decoder.GetActiveModel() != model:
                    self._SetState(self.q)
                self.local_model_steps.append(self.decoder.GetActiveModel())

    def Predict(self):
        self.Initialize()
        self._EncodeInitialState()
        self.scheme.Predict(self.strategy_data)

    def _SetState(self, q):
        """u = D(q) on the free DOFs, through the update of the scheme"""
        current = numpy.array(self.strategy_data.GetEffectiveDofSet().GetValues())
        effective_dx = self.strategy_data.GetEffectiveLinearSystem().GetVector(KratosMultiphysics.Future.DenseVectorTag.Dx)
        numpy.array(effective_dx, copy=False)[:] = self.decoder.Decode(q) - current
        self.scheme.Update(self.strategy_data)

    def _ApplyReducedIncrement(self, dq):
        """Returns the increment of the DOFs (as given by the decoder)"""
        dx = self.decoder.Decode(self.q + dq) - self.decoder.Decode(self.q)
        if self.update == "increment":
            effective_dx = self.strategy_data.GetEffectiveLinearSystem().GetVector(KratosMultiphysics.Future.DenseVectorTag.Dx)
            numpy.array(effective_dx, copy=False)[:] = dx
            self.scheme.Update(self.strategy_data)
        else:
            self._SetState(self.q + dq)
        self.q = self.q + dq
        return dx

    def _Project(self):
        """Reduced system of the assembled effective system for the projection of this strategy. Returns the reduced matrix and the
        vector whose norm measures the reduced residual (W^T b; for LSPG the gradient (A Phi)^T b of the least squares)."""
        phi = self.decoder.Jacobian(self.q)
        if self.projection == "galerkin" and not self.monotonicity_preserving:
            self.reduced_system = None
            return self.projector.Project(phi)  # C++
        # Other test bases: the assembled matrix in SciPy (as the NumpyRomProjector of rom_future_solver.py). The rows of the fixed
        # DOFs do not enter, in the trial basis nor in the test one
        effective_linear_system = self.strategy_data.GetEffectiveLinearSystem()
        lhs = scipy_conversion_tools.to_csr(effective_linear_system.GetMatrix(KratosMultiphysics.Future.SparseMatrixTag.LHS))
        free = numpy.array([not dof.IsFixed() for dof in self.strategy_data.GetEffectiveDofSet()])
        rhs = numpy.array(effective_linear_system.GetVector(KratosMultiphysics.Future.DenseVectorTag.RHS)) * free
        if self.monotonicity_preserving:
            # As MonotonicityPreserving of the legacy GlobalROMBuilderAndSolver: with S the positive entries above the diagonal and
            # their transpose, and L = diag(S 1) - S, the system becomes A + L and b - L u
            positive = scipy.sparse.triu(lhs, k=1).tocsr()
            positive.data = numpy.where(positive.data > 0.0, positive.data, 0.0)
            symmetric = positive + positive.T
            laplacian = scipy.sparse.diags(numpy.asarray(symmetric.sum(axis=1)).ravel()) - symmetric
            lhs = (lhs + laplacian).tocsr()
            rhs = (rhs - laplacian @ numpy.array(self.strategy_data.GetEffectiveDofSet().GetValues())) * free
        a_phi = (lhs @ (phi * free[:, None])) * free[:, None]  # computed before any further assembly
        if self.projection == "galerkin":
            trial = phi * free[:, None]
            self.reduced_system = (trial.T @ a_phi, trial.T @ rhs, "solve")
            return self.reduced_system[0], self.reduced_system[1]
        if self.projection == "petrov_galerkin":
            psi = self.left_basis * free[:, None]
            self.reduced_system = (psi.T @ a_phi, psi.T @ rhs, "least_squares" if psi.shape[1] != phi.shape[1] else "solve")
            return self.reduced_system[0], self.reduced_system[1]
        if self.run_hrom:
            # Second assembly, on the complementary mesh without weights, for the test basis
            if not hasattr(self, "lspg_complementary_entities"):
                self._SetLspgComplementaryMesh()
            self.scheme.SetHRomEntities(*self.lspg_complementary_entities)
            self.projector.BuildEffectiveSystem()
            complementary_lhs = scipy_conversion_tools.to_csr(effective_linear_system.GetMatrix(KratosMultiphysics.Future.SparseMatrixTag.LHS))
            test_basis = (complementary_lhs @ (phi * free[:, None])) * (free & self.lspg_selected_rows)[:, None]
            self.scheme.SetHRomEntities(*self.hrom_entities)
            self.reduced_system = (test_basis.T @ a_phi, test_basis.T @ rhs, "solve")
            return self.reduced_system[0], self.reduced_system[1]
        if self.lspg_solving_technique == "qr_decomposition":
            self.reduced_system = (a_phi, rhs, "least_squares")
        else:
            self.reduced_system = (a_phi.T @ a_phi, a_phi.T @ rhs, "solve")
        return self.reduced_system[0], a_phi.T @ rhs

    def _SolveReduced(self):
        """Increment of the reduced coordinates from the last projection"""
        if self.reduced_system is None:
            return self.projector.SolveReduced()
        reduced_lhs, reduced_rhs, technique = self.reduced_system
        if technique == "least_squares":
            return numpy.linalg.lstsq(reduced_lhs, reduced_rhs, rcond=None)[0]
        return numpy.linalg.solve(reduced_lhs, reduced_rhs)

    def _TestBasis(self):
        """The test basis on which the entity residuals of the HROM training are projected"""
        return self.left_basis if self.projection == "petrov_galerkin" else self.decoder.Jacobian(self.q)

    def _SolveWithReducedResidualCriterion(self):
        """Newton loop with the reduced residual criterion of the Future ROM (#14955): checked before each solve; the step ends without
        a further update, and the system is built once more after the last iteration to check it."""
        is_converged, initial_residual_norm, solves = False, None, 0
        for iteration in range(1, self.max_iterations + 2):
            self.model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] = iteration
            if self.update == "absolute":
                self._SetState(self.q)
            self.scheme.InitializeNonLinIteration(self.strategy_data)
            self.projector.BuildEffectiveSystem()
            _, reduced_rhs = self._Project()
            residual_norm = float(numpy.linalg.norm(reduced_rhs))
            if initial_residual_norm is None:
                initial_residual_norm = residual_norm
            epsilon = numpy.finfo(float).eps
            ratio = 0.0 if residual_norm < epsilon else residual_norm / max(initial_residual_norm, epsilon)
            is_converged = ratio <= self.residual_relative_tolerance or residual_norm / numpy.sqrt(max(len(reduced_rhs), 1)) < self.residual_absolute_tolerance
            if is_converged or iteration > self.max_iterations:
                self.scheme.FinalizeNonLinIteration(self.strategy_data)
                break
            self._ApplyReducedIncrement(self._SolveReduced())
            solves += 1
            self.scheme.FinalizeNonLinIteration(self.strategy_data)
        self.iterations.append(solves)
        self.rom_solutions.append(self.q.copy())
        return is_converged

    def _FollowSnapshot(self):
        """The step at the snapshot projected on the decoder, with no reduced system: the free DOFs move from the current state by
        the increment of the decoder closest to the snapshot (least squares on its Jacobian; Gauss-Newton for a nonlinear decoder)"""
        step = len(self.rom_solutions)
        if step >= self.snapshots.shape[1]:
            raise Exception("The snapshots have {} steps, the simulation asks for step {}.".format(self.snapshots.shape[1], step + 1))
        dof_set = self.strategy_data.GetEffectiveDofSet()
        target = numpy.asarray(self.snapshots[self.snapshot_rows, step], dtype=float)
        free = numpy.array([not dof.IsFixed() for dof in dof_set])
        self.scheme.InitializeNonLinIteration(self.strategy_data)
        for _ in range(20):
            jacobian = self.decoder.Jacobian(self.q)[free]
            dq = numpy.linalg.lstsq(jacobian, (target - numpy.array(dof_set.GetValues()))[free], rcond=None)[0]
            effective_dx = self.strategy_data.GetEffectiveLinearSystem().GetVector(KratosMultiphysics.Future.DenseVectorTag.Dx)
            numpy.array(effective_dx, copy=False)[:] = self.decoder.Decode(self.q + dq) - self.decoder.Decode(self.q)
            self.scheme.Update(self.strategy_data)
            self.q = self.q + dq
            if numpy.linalg.norm(dq) <= 1.0e-12 * max(numpy.linalg.norm(self.q), 1.0e-300):
                break
        self.scheme.FinalizeNonLinIteration(self.strategy_data)
        self.iterations.append(0)
        self.rom_solutions.append(self.q.copy())
        return True

    def SolveSolutionStep(self):
        if self.snapshots is not None:
            return self._FollowSnapshot()
        if self.convergence_criterion == "reduced_residual_criterion":
            return self._SolveWithReducedResidualCriterion()
        arguments = (self.model_part, self.dof_set, self.criterion_A, self.criterion_dx)
        iteration, is_converged = 0, False
        while not is_converged and iteration < self.max_iterations:
            iteration += 1
            self.model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] = iteration
            self.scheme.InitializeNonLinIteration(self.strategy_data)
            self.criterion.InitializeNonLinearIteration(*arguments, self.criterion_b)
            is_converged = self.criterion.PreCriteria(*arguments, self.criterion_b)

            # Full order system at the current state (Future scheme), Galerkin projection on the Jacobian of the decoder, reduced solve
            if self.update == "absolute":
                self._SetState(self.q)
            self.projector.BuildEffectiveSystem()
            reduced_lhs, _ = self._Project()
            self.max_condition_number = max(getattr(self, "max_condition_number", 0.0), numpy.linalg.cond(reduced_lhs))
            dq = self._SolveReduced()

            # Update of the decoder, applied by the scheme (free DOFs and its time derivatives)
            dx = self._ApplyReducedIncrement(dq)

            numpy.array(self.criterion_dx, copy=False)[:] = dx
            self.scheme.FinalizeNonLinIteration(self.strategy_data)
            self.criterion.FinalizeNonLinearIteration(*arguments, self.criterion_b)
            if is_converged:
                is_converged = self.criterion.PostCriteria(*arguments, self._CriterionResidual())
        self.iterations.append(iteration)
        self.rom_solutions.append(self.q.copy())
        return is_converged

    def FinalizeSolutionStep(self):
        if self.train_hrom:
            if os.environ.get("TRAIN_VIA_ENTITY") == "1":
                # Control of the study: the residual asked to the entity alone (RomProjector of #14953), without the scheme
                self.hrom_residual_snapshots.append(self.projector.ProjectEntityResiduals(self.decoder.Jacobian(self.q)))
            else:
                self.hrom_residual_snapshots.append(self.scheme.ProjectEntityResiduals(self._TestBasis(), self.strategy_data))
        self.scheme.FinalizeSolutionStep(self.strategy_data)
        self.criterion.FinalizeSolutionStep(self.model_part, self.dof_set, self.criterion_A, self.criterion_dx, self.criterion_b)

    def Check(self):
        return self.legacy_scheme.Check(self.model_part)

    def Clear(self):
        pass

    def SetEchoLevel(self, level):
        pass

    def GetScheme(self):
        return self.legacy_scheme

    # --- helpers --------------------------------------------------------------------------------------------------------------------
    def _SetHRomEntities(self):
        """HROM entities and weights of the ROM folder (as _SetHRomEntities of the FutureRomAnalysis)"""
        folder = pathlib.Path(self.rom_data_folder)
        def Weights(name, number_of_entities):
            weights = numpy.load(folder / f"HROM_{name}Weights.npy")
            if weights.size == 0:  # no entity of this kind selected
                return []
            return [float(w) for w in weights.reshape(number_of_entities, -1)[:, self.hrom_weights_set]]
        element_ids, condition_ids = numpy.load(folder / "HROM_ElementIds.npy"), numpy.load(folder / "HROM_ConditionIds.npy")
        self.hrom_entities = ([int(i) for i in element_ids], Weights("Element", len(element_ids)), [int(i) for i in condition_ids], Weights("Condition", len(condition_ids)))
        self.scheme.SetHRomEntities(*self.hrom_entities)

    def _SetLspgComplementaryMesh(self):
        """For the LSPG HROM (as the legacy LSPG builder): a row of A Phi needs every entity sharing its DOF, so the test basis is
        assembled without weights on the selected entities and their nodal neighbours (complementary mesh), and only its rows of the DOFs
        of the selected entities are used."""
        element_ids, _, condition_ids, _ = self.hrom_entities
        elements, conditions = self.model_part.Elements, self.model_part.Conditions
        element_nodes = {node.Id for i in element_ids for node in elements[i].GetGeometry()}
        condition_nodes = {node.Id for i in condition_ids for node in conditions[i].GetGeometry()}
        complementary_elements = sorted(set(element_ids) | {e.Id for e in elements if any(node.Id in element_nodes for node in e.GetGeometry())})
        complementary_conditions = sorted(set(condition_ids) | {c.Id for c in conditions if any(node.Id in element_nodes | condition_nodes for node in c.GetGeometry())})
        self.lspg_complementary_entities = (complementary_elements, [1.0] * len(complementary_elements), complementary_conditions, [1.0] * len(complementary_conditions))
        selected_nodes = element_nodes | condition_nodes
        self.lspg_selected_rows = numpy.array([dof.Id() in selected_nodes for dof in self.strategy_data.GetEffectiveDofSet()])

    def _CriterionResidual(self):
        if self.criterion_residual == "full" and self.criterion.GetActualizeRHSflag():
            rhs = self.strategy_data.GetLinearSystem().GetVector(KratosMultiphysics.Future.DenseVectorTag.RHS)
            rhs.SetValue(0.0)
            self.scheme.Build(rhs)
            numpy.array(self.criterion_b, copy=False)[:] = numpy.array(rhs, copy=False)
        return self.criterion_b


def GetStrategySolvers(solver, name=None):
    """The solvers that own a solution strategy, with their names: the solver itself, or the sub-solvers of a coupled one (any depth),
    named by their attribute in the coupled solver (e.g. "fluid_solver")."""
    sub_solvers = [(attribute, value) for attribute, value in vars(solver).items() if isinstance(value, PythonSolver)]
    if not sub_solvers:
        return [(name, solver)]
    return [leaf for attribute, sub_solver in sub_solvers for leaf in GetStrategySolvers(sub_solver, attribute)]


def UseFutureRomStrategies(analysis, rom_data_folder="rom_data", per_solver=None, **strategy_arguments):
    """Makes every solver of a (FOM) analysis stage create a FutureRomStrategy as its solution strategy. To be called before Initialize.
    strategy_arguments go to every strategy; per_solver = {solver name: arguments} overrides them for a solver (the names are the ones of
    the sub-solvers in a coupled solver, e.g. "fluid_solver"; a single solver has the name None).
    The column of the HROM weights of each solver is its position in the 'coupled_solvers' of the RomParameters.json, if given.
    Returns the solvers; their strategy is solver._GetSolutionStrategy() once the analysis is initialized."""
    with open(pathlib.Path(rom_data_folder) / "RomParameters.json") as parameter_file:
        coupled_solvers = json.load(parameter_file).get("coupled_solvers", [])
    solvers = []
    for name, solver in GetStrategySolvers(analysis._GetSolver()):
        arguments = dict(strategy_arguments, hrom_weights_set=coupled_solvers.index(name) if name in coupled_solvers else 0, name=name)
        arguments.update((per_solver or {}).get(name, {}))
        solver._CreateSolutionStrategy = lambda solver=solver, arguments=arguments: FutureRomStrategy(solver, rom_data_folder, **arguments)
        solvers.append(solver)
    # The initial state is encoded once the analysis is initialized (initial conditions applied), as the legacy enhanced ROMs do
    initialize = analysis.Initialize
    def InitializeAndEncodeInitialState():
        initialize()
        for solver in solvers:
            strategy = solver._GetSolutionStrategy()
            if strategy.initial_reduced_coordinates is not None:
                strategy.Initialize()
                strategy._EncodeInitialState()
    analysis.Initialize = InitializeAndEncodeInitialState
    return solvers


def SaveBasis(rom_data_folder, snapshots, node_ids, nodal_unknowns, number_of_modes, solver_unknowns=None):
    """Writes the ROM folder of a linear basis (RightBasisMatrix.npy, NodeIds.npy, RomParameters.json) from full order snapshots.
    snapshots: one column per snapshot, rows as node_ids times nodal_unknowns (usually increments with respect to the initial state).
    solver_unknowns: for a coupled solver, {solver name: its unknowns}: one basis of number_of_modes per solver, each in its rows
    (names as in UseFutureRomStrategies, e.g. "fluid_solver"). None: one basis for all the unknowns.
    Returns, for each basis, the first singular value left out relative to the first one."""
    folder = pathlib.Path(rom_data_folder)
    folder.mkdir(parents=True, exist_ok=True)
    snapshots = numpy.asarray(snapshots)
    basis = numpy.zeros((snapshots.shape[0], number_of_modes))
    left_out = {}
    for name, unknowns in (solver_unknowns or {None: nodal_unknowns}).items():
        rows = numpy.concatenate([numpy.arange(nodal_unknowns.index(unknown), snapshots.shape[0], len(nodal_unknowns)) for unknown in unknowns])
        modes, singular_values, _ = numpy.linalg.svd(snapshots[rows], full_matrices=False)
        basis[rows] = modes[:, :number_of_modes]
        left_out[name] = float(singular_values[number_of_modes] / singular_values[0]) if len(singular_values) > number_of_modes else 0.0
    numpy.save(folder / "RightBasisMatrix.npy", basis)
    numpy.save(folder / "NodeIds.npy", numpy.asarray(node_ids))
    rom_parameters = {"rom_format": "numpy", "projection_strategy": "galerkin",
                      "rom_settings": {"nodal_unknowns": list(nodal_unknowns), "number_of_rom_dofs": number_of_modes}}
    if solver_unknowns:
        rom_parameters["coupled_solvers"] = list(solver_unknowns)
    (folder / "RomParameters.json").write_text(json.dumps(rom_parameters, indent=4))
    return left_out


def SelectHRomEntities(solvers, rom_data_folder="rom_data", tolerance=1.0e-6, residuals=None):
    """HROM elements, conditions and weights from the training residuals, saved in the ROM folder (HROM_ElementIds.npy, ...), where
    the strategies with run_hrom=True read them. For each solver: SVD of its residuals truncated at the tolerance and empirical cubature
    (both of the RomApp, as the legacy HRomTrainingUtility), giving one column of weights per solver.
    solvers: the ones returned by UseFutureRomStrategies with train_hrom=True, after the run.
    residuals: one list of residual snapshots per solver, to train on several runs (the hrom_residual_snapshots of the strategy of
    each run, added up). None: the ones of the strategies of these solvers.
    Returns one dict per solver with the sizes."""
    folder = pathlib.Path(rom_data_folder)
    strategies = [solver._GetSolutionStrategy() for solver in solvers]
    number_of_columns = max(strategy.hrom_weights_set for strategy in strategies) + 1
    selection, report = {}, []
    for index, (solver, strategy) in enumerate(zip(solvers, strategies)):
        model_part = solver.GetComputingModelPart()
        entities = [("Element", element.Id) for element in model_part.Elements] + [("Condition", condition.Id) for condition in model_part.Conditions]
        snapshots = numpy.hstack(strategy.hrom_residual_snapshots if residuals is None else residuals[index])
        cubature = EmpiricalCubatureMethod()
        with contextlib.redirect_stdout(io.StringIO()):  # both print every iteration
            residual_modes, _, _, _ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(snapshots, tolerance)
            cubature.SetUp(residual_modes, constrain_sum_of_weights=True, number_of_conditions=model_part.NumberOfConditions())
            cubature.Run()
        selected = [entities[i] for i in numpy.squeeze(cubature.z, axis=None).reshape(-1)]
        for entity, weight in zip(selected, numpy.asarray(cubature.w).reshape(-1)):
            selection.setdefault(entity, numpy.zeros(number_of_columns))[strategy.hrom_weights_set] = weight
        report.append({"solver": strategy.name or type(solver).__name__, "residuals": snapshots.shape, "residual modes": residual_modes.shape[1],
                       "elements": sum(kind == "Element" for kind, _ in selected), "of elements": model_part.NumberOfElements(),
                       "conditions": sum(kind == "Condition" for kind, _ in selected), "of conditions": model_part.NumberOfConditions()})
    for kind in ("Element", "Condition"):
        ids = sorted(entity_id for entity_kind, entity_id in selection if entity_kind == kind)
        numpy.save(folder / f"HROM_{kind}Ids.npy", numpy.array(ids, dtype=int))
        numpy.save(folder / f"HROM_{kind}Weights.npy", numpy.array([selection[(kind, entity_id)] for entity_id in ids]).reshape(len(ids), number_of_columns))
    return report


def DescribeStrategies(solvers):
    """What stands in each slot of each solver: one dict per solver (after the analysis is initialized or run)"""
    rows = []
    for solver in solvers:
        strategy = solver._GetSolutionStrategy()
        rows.append({
            "solver": strategy.name or type(solver).__name__,
            "scheme of the solver": type(strategy.legacy_scheme).__name__,
            "Future scheme": "{} ({})".format(type(strategy.scheme).__name__, strategy.scheme_origin),
            "criterion": type(strategy.criterion).__name__ + " (solver)" if strategy.convergence_criterion == "solver" else "reduced residual (ROM)",
            "projection": strategy.projection + (" ({})".format(strategy.lspg_solving_technique) if strategy.projection == "lspg" else ""),
            "decoder": "{} {}".format(type(strategy.decoder).__name__, numpy.shape(strategy.decoder.Jacobian(strategy.q))),
            "update": strategy.update,
            "Newton solves": int(sum(strategy.iterations)),
        })
    return rows
