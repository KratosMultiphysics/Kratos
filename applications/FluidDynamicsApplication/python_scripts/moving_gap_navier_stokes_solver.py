#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   Fluid Dynamics Application
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""Experimental nonconservative Q1 Navier–Stokes moving Gap-SBM solver."""
import KratosMultiphysics as KM
from KratosMultiphysics import python_linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_stokes_solver import MovingGapStokesSolver


def create_qsvms_picard_strategy(model_part, linear_solver, settings):
    """Build the combined local system each iteration for QSVMS Picard.

    QSVMS does not implement its BDF RHS-only path; rebuilding the matrix is
    required both on the first iteration and on subsequent iterations.
    """
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    criterion = KM.DisplacementCriteria(settings['relative_tolerance'].GetDouble(),
                                        settings['absolute_tolerance'].GetDouble())
    criterion.SetEchoLevel(settings['echo_level'].GetInt())
    builder = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(
        model_part, scheme, criterion, builder, settings['maximum_iterations'].GetInt(),
        False, False, False)
    strategy.SetRebuildLevel(2)
    strategy.SetEchoLevel(settings['echo_level'].GetInt())
    return strategy


def CreateSolver(model, custom_settings):
    return MovingGapNavierStokesSolver(model, custom_settings)


class MovingGapNavierStokesSolver(MovingGapStokesSolver):
    """Nonlinear QSVMS strategy using the moving patch lifecycle."""

    @classmethod
    def GetDefaultParameters(cls):
        settings = MovingGapStokesSolver.GetDefaultParameters()
        settings['coupling_settings'].AddString('formulation', 'qsvms')
        settings.AddInt('maximum_iterations', 50)
        settings.AddDouble('relative_tolerance', 1e-8)
        settings.AddDouble('absolute_tolerance', 1e-10)
        settings.AddDouble('maximum_free_dof_residual', 1e-8)
        return settings

    def SetExactHistoryCallback(self, callback):
        if callback is not None:
            raise ValueError('Navier–Stokes moving patches require numerical history transfer.')
        super().SetExactHistoryCallback(None)

    def Initialize(self):
        coupling = self.settings['coupling_settings']
        if not coupling.Has('formulation') or coupling['formulation'].GetString() != 'qsvms':
            raise ValueError('Moving Navier–Stokes requires coupling formulation "qsvms".')
        self._converged = False
        super().Initialize()

    def InitializeSolutionStep(self):
        if self._strategy is not None:
            raise RuntimeError('The previous strategy was not cleared.')
        self._solved = False
        self._converged = False
        output = KM.Logger.GetDefaultOutput()
        severity = output.GetSeverity()
        try:
            if self.echo_level <= 1:
                output.SetSeverity(KM.Logger.Severity.WARNING)
            linear_solver = python_linear_solver_factory.ConstructSolver(
                self.settings['linear_solver_settings'].Clone())
        finally:
            output.SetSeverity(severity)
        self._strategy = create_qsvms_picard_strategy(
            self.GetComputingModelPart(), linear_solver, self.settings)
        self._strategy.Initialize()
        self._strategy.Check()

    def SolveSolutionStep(self):
        if self._strategy is None:
            raise RuntimeError('InitializeSolutionStep must precede the solve.')
        self._strategy.InitializeSolutionStep()
        self._strategy.Predict()
        converged = self._strategy.SolveSolutionStep()
        self._strategy.FinalizeSolutionStep()
        if converged:
            # QSVMS has no RHS-only time-integrated assembly. Recompute the
            # combined local systems to detect a stalled Picard iteration.
            part = self.GetComputingModelPart()
            residual = {node.Id: [0.0, 0.0, 0.0] for node in part.Nodes}
            for entity in (*part.Elements, *part.Conditions):
                lhs, rhs = KM.Matrix(), KM.Vector()
                entity.CalculateLocalSystem(lhs, rhs, part.ProcessInfo)
                for i, node in enumerate(entity.GetGeometry()):
                    for d in range(3):
                        residual[node.Id][d] += rhs[3*i+d]
            self.free_dof_residual = max((abs(residual[node.Id][d]) for node in part.Nodes
                for d, variable in enumerate((KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE))
                if not node.IsFixed(variable)), default=0.0)
            converged = self.free_dof_residual <= self.settings['maximum_free_dof_residual'].GetDouble()
        self._converged = converged
        self._solved = converged
        return converged

    def FinalizeSolutionStep(self):
        if not self._converged:
            raise RuntimeError('Cannot commit moving-patch histories before a converged solve.')
        super().FinalizeSolutionStep()
        self._converged = False
