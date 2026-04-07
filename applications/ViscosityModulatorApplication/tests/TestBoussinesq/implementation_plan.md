# Fully Implicit Coupling: Boussinesq Body Force ↔ Navier-Stokes

## Background

Currently the simulation runs two solvers sequentially each step:

1. **Navier-Stokes** (N-S) — solved with the body force `BODY_FORCE` from step n−1.
2. **Convection-Diffusion** (C-D) — solved with the velocity field from step n.

The `ApplyBoussinesqModulatorFieldProcess.ExecuteInitializeSolutionStep()` runs
**before** both solves, so the body force entering the N-S solve is always one step stale.
This is a **lagged** (explicit) coupling.

A **fully implicit** (monolithic / partitioned-implicit) coupling requires that the body
force `g·(ρ(φ)−ρ₀)/ρ₀` is consistent with the **current** concentration field φ when
the fluid solve converges.

---

## Approaches Considered

### Option A — Inject into the Element (Monolithic)

Modify `ViscosityModulatorElementShockCapturing` (or a coupled N-S element) to add the
buoyancy term directly into the element stiffness matrix so that both φ and **v**
are in the same global DOF vector. This is the most consistent approach but requires
C++ element work for both N-S and C-D simultaneously and is by far the most intrusive.

**→ Rejected** for now as it requires co-design of a monolithic element.

### Option B — Inject into `InitializeNonLinIteration` of the N-S Scheme

Override the Kratos C++ scheme and call `AssignBoussinesqForce()` inside
`InitializeNonLinIteration`. The body force would then be updated at every N-S
Newton-Raphson sub-iteration.

**→ Possible but over-coupled**: The C-D field φ is not updated within the N-S NR
iterations (they are solved at different stages), so updating `BODY_FORCE` from a
stale φ inside the N-S NR loop gains nothing.

### Option C — Partitioned Picard Loop in `SolveSolutionStep` ✅ (recommended)

Implement an **outer fixed-point (Picard) iteration** in Python, analogous to the
`ConjugateHeatTransferSolver` and `PartitionedFSIBaseSolver` patterns already present
in Kratos. Each outer iteration:

1. Assigns Boussinesq body force from `φ^k`.
2. Solves N-S → produces `v^{k+1}`.
3. Solves C-D with `v^{k+1}` → produces `φ^{k+1}`.
4. Checks `‖φ^{k+1} − φ^k‖ / ‖φ^k‖ < tol`.

This requires **no C++ changes**, only Python-level rewiring.

> [!IMPORTANT]
> This is a **partitioned** (block Gauss-Seidel) coupling. It is **not** the same as
> a monolithic solve. However, it is consistent at convergence: when the Picard loop
> converges, the Boussinesq force is exactly consistent with the current φ field.
> Convergence speed depends on the coupling strength (large buoyancy / small diffusion
> will require more Picard iterations or a relaxation/acceleration strategy).

---

## User Review Required

> [!WARNING]
> **Where does your N-S solver live?** The current `ViscosityModulatorAnalysis`
> only manages a single solver (the C-D one). For the Picard loop to work, the
> N-S solve must also be drivable from the same `SolveSolutionStep`. This plan
> assumes you already have a coupled solver that runs both N-S and C-D
> (e.g. `CoupledFluidScalarSolver` or your own variant). Please confirm.
>
> If you are running N-S and C-D as two separate analysis stages (two JSON files,
> two processes), a different coordination mechanism is needed.

> [!IMPORTANT]
> **Relaxation / acceleration**: A plain Picard iteration can be slow or divergent
> for high Rayleigh-number / convection-dominated problems. The FSI application
> solves this with an **Aitken** or **MVQN** convergence accelerator
> (`FSIApplication.convergence_accelerator_factory`). The plan below includes an
> optional relaxation coefficient. We can add Aitken acceleration in a follow-up.

---

## Proposed Changes

### Architecture

```
CoupledBoussinesqSolver             (NEW — thin orchestrator)
  ├── fluid_solver                  (existing N-S solver, e.g. Monolithic VMS)
  ├── scalar_solver                (existing stationary C-D solver)
  └── boussinesq_process            (existing BoussinesqModulatorFieldProcess)

coupling_settings:
  max_coupling_iterations : 20
  coupling_relative_tolerance : 1e-5
  relaxation_factor : 1.0
```

The new solver lives entirely in Python and orchestrates the two existing solvers.

---

### Files to Create / Modify

---

#### [NEW] `coupled_boussinesq_solver.py`

Path: `applications/ViscosityModulatorApplication/python_scripts/coupled_boussinesq_solver.py`

Full new `PythonSolver` subclass. Key method is `SolveSolutionStep`:

```python
def SolveSolutionStep(self):
    tol   = self.settings["coupling_settings"]["coupling_relative_tolerance"].GetDouble()
    omega = self.settings["coupling_settings"]["relaxation_factor"].GetDouble()
    max_it = self.settings["coupling_settings"]["max_coupling_iterations"].GetInt()

    # Store φ^k (copy current concentration to a buffer variable)
    phi_var = self._get_concentration_variable()
    self._store_old_concentration(phi_var)   # saves φ into AUX_CONCENTRATION on each node

    for it in range(1, max_it + 1):
        # 1. Update body force from current φ^k
        self.boussinesq_process.ExecuteInitializeSolutionStep()

        # 2. Solve Navier-Stokes with updated body force → v^{k+1}
        fluid_converged = self.fluid_solver.SolveSolutionStep()

        # 3. Solve Convection-Diffusion with v^{k+1} → φ^{k+1}
        scalar_converged = self.scalar_solver.SolveSolutionStep()

        # 4. Relaxation: φ^{k+1} ← ω·φ^{k+1} + (1−ω)·φ^k
        if omega < 1.0:
            self._apply_relaxation(phi_var, omega)

        # 5. Check coupling convergence: ‖φ^{k+1} − φ^k‖ / ‖φ^k‖
        rel_res = self._compute_coupling_residual(phi_var)
        KratosMultiphysics.Logger.PrintInfo(
            "::[CoupledBoussinesqSolver]::",
            f"Coupling it={it}  |Δφ|/|φ| = {rel_res:.3e}")

        if rel_res < tol:
            KratosMultiphysics.Logger.PrintInfo(
                "::[CoupledBoussinesqSolver]::",
                f"Converged in {it} iterations.")
            break

        # Save φ^{k+1} as new φ^k for next iteration
        self._store_old_concentration(phi_var)

        if it == max_it:
            KratosMultiphysics.Logger.PrintWarning(
                "::[CoupledBoussinesqSolver]::",
                "Max coupling iterations reached without convergence!")

    return fluid_converged and scalar_converged
```

Other methods delegate to the sub-solvers exactly as `CoupledFluidScalarSolver` does
(i.e. `AddVariables`, `ImportModelPart`, `PrepareModelPart`, `AddDofs`,
`InitializeSolutionStep`, `FinalizeSolutionStep`, `AdvanceInTime` all call both
sub-solvers).

The auxiliary nodal variable `AUX_CONCENTRATION` (a non-historical scalar) is used to
store `φ^k` between Picard iterations. It can be a simple `std::unordered_map<Node*, double>`
via `SetValue/GetValue` rather than a historical buffer, avoiding buffer resizing.

---

#### [MODIFY] `python_solvers_wrapper_convection_diffusion.py`

Add a new entry for your solver type:

```python
elif solver_type == "coupled_boussinesq":
    solver_module_name = "coupled_boussinesq_solver"
```

---

#### [MODIFY] `apply_boussinesq_concentration_field_process.py`

Remove the call from `ExecuteInitializeSolutionStep` (or keep it as a safety fallback
for the first step, but mark it clearly). The update is now driven explicitly by the
Picard loop in the solver.

```python
def ExecuteInitializeSolutionStep(self):
    # Body force update is now driven by CoupledBoussinesqSolver.SolveSolutionStep().
    # Calling it here would break the coupling loop. Do NOT call AssignBoussinesqForce.
    pass
```

> [!CAUTION]
> If you keep `ExecuteInitializeSolutionStep` active while also calling it inside the
> Picard loop, the body force will be overwritten with the **lagged** value at the
> start of each step, undoing the first Picard iteration. You must disable one or the
> other.

---

#### [MODIFY] `ProjectParameters.json`

Change `solver_type` to `"coupled_boussinesq"` and add coupling settings:

```json
"solver_settings": {
    "solver_type": "coupled_boussinesq",
    "fluid_solver_settings": { ... },
    "scalar_solver_settings": { ... },
    "boussinesq_process_settings": {
        "model_part_name": "FluidModelPart",
        "base_fluid_density": 1.0,
        "max_density": 2.0,
        "gravity": [0.0, -9.81, 0.0],
        "modify_pressure": false,
        "modify_density": false
    },
    "coupling_settings": {
        "max_coupling_iterations": 20,
        "coupling_relative_tolerance": 1.0e-5,
        "relaxation_factor": 1.0
    }
}
```

---

## Coupling Convergence Measure

The coupling residual is measured on the **concentration field** φ, which is the
variable that drives the body force:

```
rel_res = ‖φ^{k+1} − φ^k‖_L2 / max(‖φ^k‖_L2, 1e-10)
```

This is computed purely in Python using `FastGetSolutionStepValue`.

---

## Optional Extension: Aitken Acceleration

If plain Picard converges slowly (typical for high Ra), the relaxation factor `ω` can
be computed dynamically using the Aitken method at each coupling iteration:

```
ω^{k+1} = −ω^k · (Δφ^k · (Δφ^{k+1} − Δφ^k)) / ‖Δφ^{k+1} − Δφ^k‖²
```

This is already implemented in Kratos' `FSIApplication.ConvergenceAcceleratorFactory`
with `"solver_type": "Relaxation", "acceleration_type": "Aitken"`. We could wire it
in if needed.

---

## Verification Plan

### Functional Check
1. Run a square cavity benchmark at low Ra (e.g. Ra=10³). The Picard loop should
   converge in 1–3 iterations. Compare velocity and concentration profiles against
   the lagged solution to confirm they differ only in the expected direction.

2. Check the logger output shows decreasing `|Δφ|/|φ|` per coupling iteration.

### Correctness Check
3. At high Ra (where the lagged solution is clearly wrong), verify that the implicitly
   coupled solution gives a physically correct steady-state field.

4. Run with `relaxation_factor = 1.0` and `max_coupling_iterations = 1` — this should
   exactly replicate the current explicit (lagged) behavior, confirming backward
   compatibility.

### Convergence Rate Check
5. Plot `|Δφ|/|φ|` vs. coupling iteration number for several Ra values. If the rate
   is poor (> ~5 iterations needed), enable Aitken acceleration.

---

## Open Questions

> [!IMPORTANT]
> 1. **Do you already have a combined solver** (N-S + C-D in one `ProjectParameters.json`)?
>    If so, which solver type? (`ScalarCoupled`, your own, …?)
>
> 2. **Is your N-S solver the standard VMS Monolithic** (`navier_stokes_solver_vmsmonolithic`)
>    or a different one (e.g., fractional step)?
>
> 3. **Do you want Aitken acceleration** from the start, or start with plain Picard
>    and add it only if convergence is slow?
>
> 4. **Stationary or transient?** For a stationary problem (single time step), the
>    Picard loop is the entire solve. For a transient problem, convergence of the Picard
>    loop at each step determines the accuracy of the implicit coupling.
