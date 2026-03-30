# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import the base coupled solver that this class extends
from KratosMultiphysics.ConvectionDiffusionApplication.coupled_fluid_thermal_solver import CoupledFluidThermalSolver


def CreateSolver(model, custom_settings):
    return CoupledFluidThermalSolverBoussinesq(model, custom_settings)


class CoupledFluidThermalSolverBoussinesq(CoupledFluidThermalSolver):
    """Coupled fluid-thermal solver with fully implicit Boussinesq buoyancy coupling.

    This solver extends CoupledFluidThermalSolver by wrapping the sequential
    fluid + thermal solves inside a Picard (fixed-point) iteration loop, so that
    the Boussinesq body force computed from the concentration field φ is consistent
    with the converged state of both the Navier-Stokes and Convection-Diffusion
    problems at the end of each time step.

    Each Picard iteration k performs:
        1. Compute BODY_FORCE from current φ^k  (Boussinesq process)
        2. Solve N-S → v^{k+1}
        3. Solve C-D (convecting φ with v^{k+1}) → φ^{k+1}
        4. Optional relaxation: φ ← ω·φ^{k+1} + (1−ω)·φ^k
        5. Convergence check: ‖φ^{k+1} − φ^k‖ / ‖φ^k‖ < tol

    The Boussinesq process is owned and driven entirely by this solver. It should
    NOT appear in the 'processes' section of ProjectParameters.json — its settings
    are instead embedded under 'solver_settings["boussinesq_settings"]'.

    This design is stationary-first but fully compatible with transient runs:
    the Picard loop re-converges the coupling at every time step.

    JSON configuration example
    --------------------------
    "solver_settings": {
        "solver_type": "thermally_coupled_boussinesq",
        "domain_size": 2,
        "fluid_solver_settings":   { ... },   // unchanged from ThermallyCoupled
        "thermal_solver_settings": { ... },   // unchanged from ThermallyCoupled
        "boussinesq_settings": {
            "base_fluid_density": 1.0,
            "max_density": 2.0,
            "gravity": [0.0, -9.81, 0.0],
            "modify_pressure": false,
            "modify_density": false,
            "r0": [0.0, 0.0, 0.0]
        },
        "coupling_settings": {
            "max_coupling_iterations": 20,
            "coupling_relative_tolerance": 1.0e-5,
            "relaxation_factor": 1.0
        }
    }
    """

    # ------------------------------------------------------------------
    # Class-level defaults
    # ------------------------------------------------------------------

    @classmethod
    def GetDefaultParameters(cls):
        """Return default solver parameters, extending the base class defaults."""

        # New parameters introduced by this solver
        my_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "thermally_coupled_boussinesq",
            "boussinesq_settings": {
                "base_fluid_density" : 1.0,
                "max_density"        : 2.0,
                "gravity"            : [0.0, 0.0, 0.0],
                "modify_pressure"    : false,
                "modify_density"     : false,
                "r0"                 : [0.0, 0.0, 0.0]
            },
            "coupling_settings": {
                "max_coupling_iterations"    : 20,
                "coupling_relative_tolerance": 1.0e-5,
                "relaxation_factor"          : 1.0
            }
        }""")

        # Merge with the base class defaults (fluid/thermal solver settings, etc.)
        my_defaults.AddMissingParameters(super().GetDefaultParameters())
        return my_defaults

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, model, custom_settings):
        # Let the base class create the fluid and thermal sub-solvers
        super().__init__(model, custom_settings)

        # Cache coupling parameters for quick access at solve time
        coupling = self.settings["coupling_settings"]
        self._max_coupling_it = coupling["max_coupling_iterations"].GetInt()
        self._coupling_tol    = coupling["coupling_relative_tolerance"].GetDouble()
        self._omega           = coupling["relaxation_factor"].GetDouble()

        KratosMultiphysics.Logger.PrintInfo(
            "::[CoupledFluidThermalSolverBoussinesq]::",
            "Construction finished. "
            "Max Picard iterations: {0}, "
            "tolerance: {1:.1e}, "
            "relaxation ω: {2:.2f}".format(
                self._max_coupling_it, self._coupling_tol, self._omega))

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def Initialize(self):
        """Initialize both sub-solvers and then create the Boussinesq process.

        The process is created AFTER the base Initialize() so that:
        - The thermal model part is fully prepared and its ProcessInfo contains
          CONVECTION_DIFFUSION_SETTINGS (needed to locate the unknown variable).
        - The variable lists of both model parts are already merged, so
          BODY_FORCE is accessible on the thermal model part's nodes.
        """
        # Base class initializes fluid and thermal strategies
        super().Initialize()

        # Build the Boussinesq process, pointing it at the thermal model part.
        # The thermal model part has CONVECTION_DIFFUSION_SETTINGS in its ProcessInfo
        # and shares nodes with the fluid model part (via ConnectivityPreserveModeler),
        # so modifying BODY_FORCE here affects the fluid solve directly.
        boussinesq_settings = self.settings["boussinesq_settings"]
        self._boussinesq_process = ConvectionDiffusionApplication.BoussinesqConcentrationFieldProcess(
            self.thermal_solver.main_model_part,
            boussinesq_settings)
        self._boussinesq_process.ExecuteInitialize()

        KratosMultiphysics.Logger.PrintInfo(
            "::[CoupledFluidThermalSolverBoussinesq]::",
            "Boussinesq process initialized.")

    # ------------------------------------------------------------------
    # Per-step solve — Picard coupling loop
    # ------------------------------------------------------------------

    def SolveSolutionStep(self):
        """Solve both physics with fully implicit Boussinesq coupling.

        Runs a Picard (block Gauss-Seidel) loop that alternately updates
        the Boussinesq body force, solves N-S, and solves C-D until the
        concentration field converges between iterations.

        Returns
        -------
        bool
            True if both sub-solvers reported convergence in the last
            Picard iteration.
        """
        # Retrieve the concentration (unknown) variable from the C-D settings
        conv_diff_settings = self.thermal_solver.main_model_part.ProcessInfo[
            KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS]
        phi_var = conv_diff_settings.GetUnknownVariable()

        # φ^0 = current concentration on nodes (from previous step / initial condition)
        phi_old = self._store_concentration(phi_var)

        fluid_converged   = True
        thermal_converged = True
        is_coupling_converged = False
        rel_res = 0.0

        for it in range(1, self._max_coupling_it + 1):

            # --- Step 1: update body force from current φ -------------------
            # Uses φ on the nodes (= φ^{k-1} at the start of iteration k).
            self._boussinesq_process.Execute()

            # --- Step 2: solve Navier-Stokes → v^k -------------------------
            fluid_converged = self.fluid_solver.SolveSolutionStep()

            # --- Step 3: solve Convection-Diffusion → φ^k ------------------
            # The C-D solver reads VELOCITY from the now-updated fluid model part.
            thermal_converged = self.thermal_solver.SolveSolutionStep()

            # --- Step 4: relaxation (only if ω < 1) -------------------------
            phi_new = self._store_concentration(phi_var)
            if self._omega < 1.0:
                self._apply_relaxation(phi_var, phi_old, phi_new, self._omega)
                phi_new = self._store_concentration(phi_var)

            # --- Step 5: coupling convergence check -------------------------
            rel_res = self._compute_residual(phi_old, phi_new)
            KratosMultiphysics.Logger.PrintInfo(
                "::[CoupledFluidThermalSolverBoussinesq]::",
                "Picard it={0:3d} | |Δφ|/|φ| = {1:.3e}".format(it, rel_res))

            # φ^{k+1} becomes φ^k for the next iteration
            phi_old = phi_new

            if rel_res <= self._coupling_tol:
                is_coupling_converged = True
                KratosMultiphysics.Logger.PrintInfo(
                    "::[CoupledFluidThermalSolverBoussinesq]::",
                    "Coupling converged in {0} Picard iteration(s).".format(it))
                break

        if not is_coupling_converged:
            KratosMultiphysics.Logger.PrintWarning(
                "::[CoupledFluidThermalSolverBoussinesq]::",
                "Picard loop reached max iterations ({0}) without convergence. "
                "Last |Δφ|/|φ| = {1:.3e}".format(self._max_coupling_it, rel_res))

        return fluid_converged and thermal_converged

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _store_concentration(self, phi_var):
        """Return a Python list containing the current nodal φ values (step 0).

        Iterates over the thermal model part nodes. Because the thermal and
        fluid model parts share the same node objects, this represents the
        concentration on the full mesh.
        """
        return [node.GetSolutionStepValue(phi_var)
                for node in self.thermal_solver.main_model_part.Nodes]

    def _compute_residual(self, phi_old, phi_new):
        """Compute the relative L2 norm of the coupling residual.

        residual = ‖φ_new − φ_old‖_L2 / max(‖φ_old‖_L2, 1e-10)

        Parameters
        ----------
        phi_old : list[float]
            Nodal φ values from the previous Picard iteration.
        phi_new : list[float]
            Nodal φ values from the current Picard iteration.

        Returns
        -------
        float
            Relative coupling residual.
        """
        diff_norm_sq = sum((a - b) ** 2 for a, b in zip(phi_new, phi_old))
        old_norm_sq  = sum(a ** 2 for a in phi_old)
        return diff_norm_sq ** 0.5 / max(old_norm_sq ** 0.5, 1.0e-10)

    def _apply_relaxation(self, phi_var, phi_old, phi_new, omega):
        """Apply under-relaxation to the concentration field on all nodes.

        Sets each node's φ ← ω·φ_new + (1−ω)·φ_old.

        Parameters
        ----------
        phi_var : Variable[double]
            The nodal scalar unknown variable (e.g. TEMPERATURE / CONCENTRATION).
        phi_old : list[float]
            Nodal values before this Picard iteration.
        phi_new : list[float]
            Nodal values after the C-D solve of this Picard iteration.
        omega : float
            Relaxation factor in (0, 1]. omega = 1.0 → no relaxation.
        """
        one_minus_omega = 1.0 - omega
        for node, f_old, f_new in zip(
                self.thermal_solver.main_model_part.Nodes, phi_old, phi_new):
            node.SetSolutionStepValue(
                phi_var, 0, omega * f_new + one_minus_omega * f_old)
