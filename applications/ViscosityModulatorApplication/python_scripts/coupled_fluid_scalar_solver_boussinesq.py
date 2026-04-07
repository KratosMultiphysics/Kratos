# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ViscosityModulatorApplication as ViscosityModulatorApplication

# Import the base coupled solver that this class extends
from KratosMultiphysics.ViscosityModulatorApplication.coupled_fluid_scalar_solver import CoupledFluidScalarSolver
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory


def CreateSolver(model, custom_settings):
    return CoupledFluidScalarSolverBoussinesq(model, custom_settings)


class CoupledFluidScalarSolverBoussinesq(CoupledFluidScalarSolver):
    """Coupled fluid-thermal solver with fully implicit Boussinesq buoyancy coupling.

    This solver extends CoupledFluidScalarSolver by wrapping the sequential
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
        "solver_type": "scalar_coupled_boussinesq",
        "domain_size": 2,
        "fluid_solver_settings":   { ... },   // unchanged from ScalarCoupled
        "scalar_solver_settings": { ... },   // unchanged from ScalarCoupled
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
            "solver_type": "scalar_coupled_boussinesq",
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
                "relaxation_factor"          : 1.0,
                "coupling_strategy"          : "Picard",
                "picard_iterations_before_newton": 0,
                "quasi_newton_settings"      : {
                    "solver_type" : "MVQN_recursive",
                    "w_0"         : 0.825
                }
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
        self._coupling_strategy = coupling["coupling_strategy"].GetString()
        
        if self._coupling_strategy == "QuasiNewton":
            self._qn_settings = coupling["quasi_newton_settings"]
            self._picard_iters_first = coupling["picard_iterations_before_newton"].GetInt()

        KratosMultiphysics.Logger.PrintInfo(
            "::[CoupledFluidScalarSolverBoussinesq]::",
            "Construction finished. "
            "Strategy: {0}, "
            "Max iterations: {1}, "
            "tolerance: {2:.1e}, "
            "relaxation ω: {3:.2f}".format(
                self._coupling_strategy, self._max_coupling_it, self._coupling_tol, self._omega))

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def Initialize(self):
        """Initialize both sub-solvers and then create the Boussinesq process.

        The process is created AFTER the base Initialize() so that:
                - The thermal model part is fully prepared and its ProcessInfo contains
                    VISCOSITY_MODULATOR_SETTINGS (needed to locate the unknown variable).
        - The variable lists of both model parts are already merged, so
          BODY_FORCE is accessible on the thermal model part's nodes.
        """
        # Base class initializes fluid and thermal strategies
        super().Initialize()

        # Build the Boussinesq process, pointing it at the thermal model part.
        # The thermal model part has VISCOSITY_MODULATOR_SETTINGS in its ProcessInfo
        # and shares nodes with the fluid model part (via ConnectivityPreserveModeler),
        # so modifying BODY_FORCE here affects the fluid solve directly.
        boussinesq_settings = self.settings["boussinesq_settings"]
        self._boussinesq_process = ViscosityModulatorApplication.BoussinesqModulatorFieldProcess(
            self.scalar_solver.main_model_part,
            boussinesq_settings)
        self._boussinesq_process.ExecuteInitialize()

        # Create Quasi-Newton convergence accelerator if requested
        if self._coupling_strategy == "QuasiNewton":
            self._convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self._qn_settings)
            self._convergence_accelerator.Initialize()

        KratosMultiphysics.Logger.PrintInfo(
            "::[CoupledFluidScalarSolverBoussinesq]::",
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
        conv_diff_settings = self.scalar_solver.main_model_part.ProcessInfo[
            ViscosityModulatorApplication.VISCOSITY_MODULATOR_SETTINGS]
        phi_var = conv_diff_settings.GetUnknownVariable()

        # φ^0 = current concentration on nodes (from previous step / initial condition)
        # Store securely inside non-historical variable AUX_SCALAR
        aux_phi_var = ViscosityModulatorApplication.AUX_SCALAR
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
            phi_var, aux_phi_var, self.scalar_solver.main_model_part, self.scalar_solver.main_model_part, 0)

        fluid_converged   = True
        scalar_converged = True
        is_coupling_converged = False
        rel_res = 0.0

        if self._coupling_strategy == "QuasiNewton":
            self._convergence_accelerator.InitializeSolutionStep()

        self._boussinesq_process.ExecuteInitializeSolutionStep()
        for it in range(1, self._max_coupling_it + 1):

            # --- Initialize QN iteration if past the Picard phase -----------
            if self._coupling_strategy == "QuasiNewton" and it > self._picard_iters_first:
                self._convergence_accelerator.InitializeNonLinearIteration()

            # --- Step 1: update body force from current φ -------------------
            # Uses φ on the nodes (= φ^{k-1} at the start of iteration k).
            self._boussinesq_process.Execute()

            # --- Step 2: solve Navier-Stokes → v^k -------------------------
            fluid_converged = self.fluid_solver.SolveSolutionStep()

            # --- Step 3: solve Convection-Diffusion → φ^k ------------------
            # The C-D solver reads VELOCITY from the now-updated fluid model part.
            scalar_converged = self.scalar_solver.SolveSolutionStep()

            # --- Step 4/5: Relaxation/Quasi-Newton and residual --------
            # Compute relative residual BEFORE relaxation/QN so we see the raw prediction delta
            rel_res = ViscosityModulatorApplication.ViscosityModulatorCouplingUtilities.ComputeRelativeResidual(
                self.scalar_solver.main_model_part, phi_var, aux_phi_var)

            KratosMultiphysics.Logger.PrintInfo(
                "::[CoupledFluidScalarSolverBoussinesq]::",
                "Coupling it={0:3d} | |Δφ|/|φ| = {1:.3e}".format(it, rel_res))

            if rel_res <= self._coupling_tol:
                is_coupling_converged = True
                KratosMultiphysics.Logger.PrintInfo(
                    "::[CoupledFluidScalarSolverBoussinesq]::",
                    "Coupling converged in {0} iteration(s).".format(it))
                break

            if self._coupling_strategy == "QuasiNewton" and it > self._picard_iters_first:
                num_nodes = self.scalar_solver.main_model_part.NumberOfNodes()
                r_vec = KratosMultiphysics.Vector(num_nodes)
                x_vec = KratosMultiphysics.Vector(num_nodes)
                
                # Build Kratos Vectors for the accelerator using C++ utility
                ViscosityModulatorApplication.ViscosityModulatorCouplingUtilities.ComputeQuasiNewtonUpdateVectors(
                    self.scalar_solver.main_model_part, phi_var, aux_phi_var, x_vec, r_vec)
                
                # Apply Newton-Raphson update
                self._convergence_accelerator.UpdateSolution(r_vec, x_vec)
                
                # Write back converged/relaxed values using C++ utility
                ViscosityModulatorApplication.ViscosityModulatorCouplingUtilities.UpdateConvergenceVariables(
                    self.scalar_solver.main_model_part, phi_var, aux_phi_var, x_vec)
                
                self._convergence_accelerator.FinalizeNonLinearIteration()
                
            else:
                # Classic Picard with fixed relaxation
                if self._omega < 1.0:
                    ViscosityModulatorApplication.ViscosityModulatorCouplingUtilities.ApplyRelaxation(
                        self.scalar_solver.main_model_part, phi_var, aux_phi_var, self._omega)

                # φ^{k+1} becomes φ^k for the next iteration
                KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                    phi_var, aux_phi_var, self.scalar_solver.main_model_part, self.scalar_solver.main_model_part, 0)

        self._boussinesq_process.ExecuteFinalizeSolutionStep()

        if not is_coupling_converged:
            KratosMultiphysics.Logger.PrintWarning(
                "::[CoupledFluidScalarSolverBoussinesq]::",
                "Coupling loop reached max iterations ({0}) without convergence. "
                "Last |Δφ|/|φ| = {1:.3e}".format(self._max_coupling_it, rel_res))

        return fluid_converged and scalar_converged


