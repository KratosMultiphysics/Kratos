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

"""Transient linear Stokes solver on a rebuilt, moving Gap-SBM model part.

The caller prepares the patches, initial fields, coupling and ``Computing``
model part before Initialize. Boundary conditions and forcing may be assigned
between AdvanceInTime and SolveSolutionStep, as in Kratos analysis stages.
"""
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.LinearSolversApplication
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics import python_linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_lifecycle import MovingGapLifecycle


def CreateSolver(model, custom_settings):
    return MovingGapStokesSolver(model, custom_settings)


class MovingGapStokesSolver(PythonSolver):
    """Rebuild and solve prepared, rigidly translating Q1 patch pairs.

    The time step is fixed; this solver currently supports linear Stokes only.
    It recreates its strategy after every rebuild because the old strategy owns
    references to elements and conditions removed during the rebuild.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return KM.Parameters('''{
            "echo_level": 0,
            "time_step": 0.01,
            "translation_velocity": [0.0, 0.0],
            "motion": {"type": "constant", "amplitude": 0.0, "frequency": 0.0},
            "cylinder": {"enabled": false, "patch_index": 0, "center": [0.5, 0.5],
                         "radius": 0.1, "clearance": 0.0, "penalty_coefficient": 0.0},
            "patches": [],
            "bdf1_startup": true,
            "geometry_model_part_name": "Interfaces",
            "computing_model_part_name": "Computing",
            "tolerance": 1e-10,
            "history_diagnostics": {"enabled": false, "output_directory": "output/history_diagnostics"},
            "coupling_settings": {
                "background_model_part_name": "Fluid.Patch_0",
                "inner_model_part_name": "Fluid.Patch_1",
                "interface_model_part_name": "Fluid.Coupling",
                "properties_id": 100,
                "tolerance": 1e-10,
                "allow_transient": true
            },
            "linear_solver_settings": {"solver_type": "LinearSolversApplication.sparse_lu"}
        }''')

    def __init__(self, model, settings):
        super().__init__(model, settings)
        self._lifecycle = None
        self._strategy = None
        self._history_callback = None
        self.motion_diagnostics = None
        self._solved = False

    def SetExactHistoryCallback(self, callback):
        """Set an optional verification callback(node, time) for missing histories."""
        if callback is not None and not callable(callback):
            raise TypeError('History callback must be callable or None.')
        self._history_callback = callback

    def Initialize(self):
        if self._lifecycle is not None:
            raise RuntimeError('The solver is already initialized.')
        patches = self.settings['patches']
        if patches.size():
            coupling_settings = [patches[i]['coupling_settings'] for i in range(patches.size())]
            velocities = [list(patches[i]['translation_velocity'].GetVector()) for i in range(patches.size())]
            motion_settings = [patches[i]['motion'] if patches[i].Has('motion') else None
                               for i in range(patches.size())]
        else:
            coupling_settings = self.settings['coupling_settings']
            velocities = list(self.settings['translation_velocity'].GetVector())
            motion_settings = [self.settings['motion']]
        motions = []
        for motion in motion_settings:
            if motion is None or motion['type'].GetString() == 'constant':
                motions.append(None)
            elif motion['type'].GetString() == 'sinusoidal':
                motions.append({'type': 'sinusoidal', 'amplitude': motion['amplitude'].GetDouble(),
                                'frequency': motion['frequency'].GetDouble()})
            else:
                raise ValueError('Motion type must be constant or sinusoidal.')
        cylinder = None
        cylinder_settings = self.settings['cylinder']
        if cylinder_settings['enabled'].GetBool():
            from KratosMultiphysics.FluidDynamicsApplication.cylinder_sbm_geometry import CylinderSbmGeometry
            from KratosMultiphysics.FluidDynamicsApplication.rebuild_fluid_patch_gap_process import RebuildFluidPatchGapProcess
            pairs = coupling_settings if isinstance(coupling_settings, list) else [coupling_settings]
            index = cylinder_settings['patch_index'].GetInt()
            if index < 0 or index >= len(pairs):
                raise ValueError('Cylinder patch index is out of range.')
            inner = self.model.GetModelPart(pairs[index]['inner_model_part_name'].GetString())
            cylinder = CylinderSbmGeometry(inner, list(cylinder_settings['center'].GetVector()),
                cylinder_settings['radius'].GetDouble(), cylinder_settings['clearance'].GetDouble(),
                cylinder_settings['penalty_coefficient'].GetDouble())
            rebuild = RebuildFluidPatchGapProcess(self.model, pairs,
                geometry_name=self.settings['geometry_model_part_name'].GetString(),
                computing_name=self.settings['computing_model_part_name'].GetString(),
                tolerance=self.settings['tolerance'].GetDouble())
            rebuild.Execute()
        self._lifecycle = MovingGapLifecycle(
            self.model, coupling_settings, velocities,
            self.settings['time_step'].GetDouble(),
            tolerance=self.settings['tolerance'].GetDouble(),
            bdf1_startup=self.settings['bdf1_startup'].GetBool(),
            geometry_name=self.settings['geometry_model_part_name'].GetString(),
            computing_name=self.settings['computing_model_part_name'].GetString(),
            motions=motions,
            cylinder=cylinder,
            history_diagnostics=(self.settings['history_diagnostics']['output_directory'].GetString()
                                 if self.settings['history_diagnostics']['enabled'].GetBool() else None))

    def AdvanceInTime(self, current_time):
        if self._lifecycle is None:
            raise RuntimeError('Initialize the solver before advancing time.')
        self.Clear()
        time = current_time + self.settings['time_step'].GetDouble()
        _, self.motion_diagnostics = self._lifecycle.prepare_solution_step(
            time, exact_history=self._history_callback)
        return time

    def InitializeSolutionStep(self):
        if self._strategy is not None:
            raise RuntimeError('The previous strategy was not cleared.')
        self._solved = False
        # The factory logs independently of the strategy's echo level.
        log_output = KM.Logger.GetDefaultOutput()
        previous_severity = log_output.GetSeverity()
        try:
            if self.echo_level <= 1:
                log_output.SetSeverity(KM.Logger.Severity.WARNING)
            linear_solver = python_linear_solver_factory.ConstructSolver(
                self.settings['linear_solver_settings'].Clone())
        finally:
            log_output.SetSeverity(previous_severity)
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        builder = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
        self._strategy = KM.ResidualBasedLinearStrategy(
            self.GetComputingModelPart(), scheme, builder, False, False, False, False)
        self._strategy.SetEchoLevel(self.echo_level)
        self._strategy.Initialize()
        self._strategy.Check()

    def SolveSolutionStep(self):
        if self._strategy is None:
            raise RuntimeError('InitializeSolutionStep must precede the solve.')
        # Solve performs the static strategy's complete per-step lifecycle.
        self._strategy.Solve()
        self._solved = True
        return True

    def FinalizeSolutionStep(self):
        if self._strategy is None or not self._solved:
            raise RuntimeError('A successful solve is required before finalizing.')
        self.Clear()
        self._lifecycle.finalize_solution_step(
            self.GetComputingModelPart().ProcessInfo[KM.TIME])
        self._solved = False

    def Clear(self):
        if self._strategy is not None:
            self._strategy.Clear()
            self._strategy = None

    def GetComputingModelPart(self):
        if self._lifecycle is None:
            raise RuntimeError('Initialize the solver first.')
        return self._lifecycle.computing
