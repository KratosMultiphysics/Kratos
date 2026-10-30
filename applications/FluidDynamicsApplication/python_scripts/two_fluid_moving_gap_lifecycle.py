#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   License: BSD License (kratos/license.txt)
#
"""Experimental fixed-grid two-fluid level-set and QSVMS lifecycle."""
import json
import math

import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_navier_stokes_solver import create_qsvms_picard_strategy
from KratosMultiphysics.FluidDynamicsApplication.two_fluid_level_set_geometry import (
    create_complete_phase_view, extract_height_contour)
from KratosMultiphysics.FluidDynamicsApplication.two_fluid_phase_history import (
    PhaseHistory, make_current_phase_velocity_sampler, transfer_physical_velocity_histories)


class TwoFluidMovingGapLifecycle:
    """Rebuild two independent Q1 phase views and solve with BDF1/BDF2.

    The two input patches must contain complete, independent Cartesian grids
    with preassigned materials, DOFs and wall constraints. The level-set
    transport must already be initialized at the fluid model's current time.
    This is an experimental nonconservative formulation: no advective
    connector flux or moving-domain geometric conservation term is supplied.
    A failed operation leaves the model unusable until restored by the caller.
    """

    def __init__(self, model, water, air, transport, dt, coupling_properties_id, *, tolerance=1e-10):
        if not math.isfinite(dt) or dt <= 0.0 or not math.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("Time step and tolerance must be finite and positive")
        if water.GetRootModelPart() != air.GetRootModelPart() or water == air:
            raise ValueError("Water and air must be independent patches of the same root")
        self.model, self.water, self.air = model, water, air
        self.root = water.GetRootModelPart()
        if self.root.GetBufferSize() < 3:
            raise ValueError("Two-fluid BDF requires a three-step buffer")
        self.transport, self.dt, self.tolerance = transport, dt, tolerance
        self.properties_id = coupling_properties_id
        if not self.root.HasProperties(coupling_properties_id):
            raise ValueError("Coupling properties must exist before initializing the lifecycle")
        if not math.isclose(self.root.ProcessInfo[KM.TIME], transport.part.ProcessInfo[KM.TIME],
                            rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("Fluid and level-set initial times must agree")
        self._names = ("TwoFluidGapCoupling", "TwoFluidGapGeometry", "TwoFluidGapComputing")
        if any(self.root.HasSubModelPart(name) for name in self._names):
            raise ValueError("Two-fluid lifecycle output model parts already exist")
        if water.HasSubModelPart("TwoFluidActive") or air.HasSubModelPart("TwoFluidActive"):
            raise ValueError("TwoFluidActive views already exist")
        self.histories = (PhaseHistory("water"), PhaseHistory("air"))
        self._last_time = self.root.ProcessInfo[KM.TIME]
        self._step = 0
        self._failed = True
        self._pending = False
        self.root.ProcessInfo[KM.DELTA_TIME] = dt
        self.root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([1.0 / dt, -1.0 / dt, 0.0])
        contour = extract_height_contour(transport.part, tolerance)
        self._build_geometry(contour)
        self._capture(self._last_time)
        self._failed = False

    def _remove_generated(self):
        coupling, geometry, computing = self._names
        if self.root.HasSubModelPart(computing):
            self.root.RemoveSubModelPart(computing)
        if self.root.HasSubModelPart(coupling):
            part = self.root.GetSubModelPart(coupling)
            for element in list(part.Elements):
                self.root.RemoveElementFromAllLevels(element.Id)
            for condition in list(part.Conditions):
                self.root.RemoveConditionFromAllLevels(condition.Id)
            self.root.RemoveSubModelPart(coupling)
        if self.root.HasSubModelPart(geometry):
            part = self.root.GetSubModelPart(geometry)
            for item in list(part.Geometries):
                self.root.RemoveGeometryFromAllLevels(item.Id)
            self.root.RemoveSubModelPart(geometry)
        for patch in (self.water, self.air):
            if patch.HasSubModelPart("TwoFluidActive"):
                patch.RemoveSubModelPart("TwoFluidActive")

    def _build_geometry(self, contour):
        self._remove_generated()
        self.phases = tuple(create_complete_phase_view(patch, "TwoFluidActive", contour, i == 0, self.tolerance)
                            for i, patch in enumerate((self.water, self.air)))
        for phase in self.phases:
            for element in phase.Elements:
                element.Set(KM.ACTIVE, True)
        matrix = KM.Matrix(len(contour), 2)
        for i, (x, y) in enumerate(contour):
            matrix[i, 0], matrix[i, 1] = x, y
        geometry = self.root.CreateSubModelPart(self._names[1])
        self.geometry = tuple(geometry.CreateSubModelPart(name) for name in ("Water", "Air"))
        for phase, view in enumerate(self.phases):
            KM.FemMultipatchGapUtilities.CreateLevelSetPhaseGeometry(
                view, self.geometry[phase], matrix, phase == 0, self.tolerance)
        settings = KM.Parameters(json.dumps({
            "background_model_part_name": self.phases[0].FullName(),
            "inner_model_part_name": self.phases[1].FullName(),
            "interface_model_part_name": self.root.FullName() + "." + self._names[0],
            "properties_id": self.properties_id,
            "formulation": "two_phase_level_set_qsvms",
            "allow_transient": True,
            "tolerance": self.tolerance,
            "level_set_contour": contour,
        }))
        CFD.CreateFluidPatchGapCouplingConditionsProcess(self.model, settings).Execute()
        coupling = self.root.GetSubModelPart(self._names[0])
        self.computing = self.root.CreateSubModelPart(self._names[2])
        self.computing.AddElements([element.Id for view in self.phases for element in view.Elements])
        self.computing.AddElements([element.Id for element in coupling.Elements])
        self.computing.AddConditions([condition.Id for condition in coupling.Conditions])
        nodes = {node.Id for entity in (*self.computing.Elements, *self.computing.Conditions)
                 for node in entity.GetGeometry()}
        self.computing.AddNodes(sorted(nodes))
        if sum(node.IsFixed(KM.PRESSURE) for node in self.computing.Nodes) != 1:
            raise ValueError("The active two-fluid system requires one global pressure reference")

    def _capture(self, time):
        for history, phase, geometry in zip(self.histories, self.phases, self.geometry):
            history.capture(phase, geometry, time)

    def prepare_solution_step(self, time):
        """Advance the level set and rebuild. Clear any old strategy first."""
        if self._failed:
            raise RuntimeError("A previous operation failed; restore a checkpoint before retrying")
        if self._pending:
            raise RuntimeError("Solve and finalize the pending step first")
        if not math.isfinite(time) or not math.isclose(time - self._last_time, self.dt, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("Two-fluid moving gap requires successive fixed time steps")
        self._failed = True
        sampler = make_current_phase_velocity_sampler(*self.histories, self.tolerance)
        contour = self.transport.advance(time, sampler)
        self.root.CloneTimeStep(time)
        self._step += 1
        for history in self.histories:
            history.advance()
        self.root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector(
            [1.0 / self.dt, -1.0 / self.dt, 0.0] if self._step == 1 else
            [1.5 / self.dt, -2.0 / self.dt, 0.5 / self.dt])
        self._build_geometry(contour)
        required = (1,) if self._step == 1 else (1, 2)
        counts = transfer_physical_velocity_histories(
            self.histories, [list(phase.Nodes) for phase in self.phases], time, self.dt, required, self.tolerance)
        for phase in self.phases:
            for node in phase.Nodes:
                node.SetSolutionStepValue(KM.VELOCITY, node.GetSolutionStepValue(KM.VELOCITY, 1))
        self._pending, self._failed = True, False
        return self.computing, counts

    def solve_solution_step(self, linear_solver, settings, before_solve=None):
        """Solve and commit only after convergence and a free-residual check."""
        if self._failed or not self._pending:
            raise RuntimeError("Prepare a valid solution step before solving")
        self._failed = True
        if before_solve is not None:
            before_solve(self.computing, self.root.ProcessInfo[KM.TIME])
        strategy = create_qsvms_picard_strategy(self.computing, linear_solver, settings)
        try:
            strategy.Initialize()
            strategy.Check()
            strategy.InitializeSolutionStep()
            strategy.Predict()
            converged = strategy.SolveSolutionStep()
            strategy.FinalizeSolutionStep()
            if not converged:
                raise RuntimeError("Two-fluid QSVMS nonlinear strategy did not converge")
            residual = {node.Id: [0.0, 0.0, 0.0] for node in self.computing.Nodes}
            for entity in (*self.computing.Elements, *self.computing.Conditions):
                lhs, rhs = KM.Matrix(), KM.Vector()
                entity.CalculateLocalSystem(lhs, rhs, self.computing.ProcessInfo)
                for i, node in enumerate(entity.GetGeometry()):
                    for d in range(3):
                        residual[node.Id][d] += rhs[3 * i + d]
            free_residual = max((abs(residual[node.Id][d]) for node in self.computing.Nodes
                                 for d, variable in enumerate((KM.VELOCITY_X, KM.VELOCITY_Y, KM.PRESSURE))
                                 if not node.IsFixed(variable)), default=0.0)
            if not math.isfinite(free_residual) or free_residual > settings["maximum_free_dof_residual"].GetDouble():
                raise RuntimeError(f"Two-fluid QSVMS free residual is {free_residual}")
            self._capture(self.root.ProcessInfo[KM.TIME])
        finally:
            strategy.Clear()
        self._last_time = self.root.ProcessInfo[KM.TIME]
        self._pending, self._failed = False, False
        return free_residual
