# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import the base coupled solver that this class extends
from KratosMultiphysics.ConvectionDiffusionApplication.coupled_fluid_thermal_solver import CoupledFluidThermalSolver
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory


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
        },
        "sub_model_part_conditions": [
            {
                "sub_model_part_name": "Walls",
                "condition_name": "ThermalFace"
            },
            {
                "sub_model_part_name": "Outlet",
                "condition_name": "ConsistentFluxBoundaryCondition"
            }
        ]
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
                "relaxation_factor"          : 1.0,
                "coupling_strategy"          : "Picard",
                "picard_iterations_before_newton": 0,
                "quasi_newton_settings"      : {
                    "solver_type" : "MVQN_recursive",
                    "w_0"         : 0.825
                }
            },
            "sub_model_part_conditions": []
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

        self._ValidateCustomConditionAssignments()

        KratosMultiphysics.Logger.PrintInfo(
            "::[CoupledFluidThermalSolverBoussinesq]::",
            "Construction finished. "
            "Strategy: {0}, "
            "Max iterations: {1}, "
            "tolerance: {2:.1e}, "
            "relaxation ω: {3:.2f}".format(
                self._coupling_strategy, self._max_coupling_it, self._coupling_tol, self._omega))

    def PrepareModelPart(self):
        if not self._UsesCustomConditionAssignments():
            super().PrepareModelPart()
            return

        self.fluid_solver.PrepareModelPart()
        self._PrepareThermalModelPartWithCustomConditions()

    def _ValidateCustomConditionAssignments(self):
        default_assignment_settings = KratosMultiphysics.Parameters("""{
            "sub_model_part_name": "",
            "condition_name": ""
        }""")

        custom_assignments = self.settings["sub_model_part_conditions"]
        assigned_sub_model_parts = set()

        for i in range(custom_assignments.size()):
            assignment = custom_assignments[i]
            assignment.ValidateAndAssignDefaults(default_assignment_settings)

            sub_model_part_name = assignment["sub_model_part_name"].GetString()
            condition_name = assignment["condition_name"].GetString()

            if sub_model_part_name == "":
                raise Exception(
                    "Each entry in 'sub_model_part_conditions' requires a non-empty 'sub_model_part_name'.")

            if condition_name == "":
                raise Exception(
                    "Each entry in 'sub_model_part_conditions' requires a non-empty 'condition_name'.")

            if sub_model_part_name in assigned_sub_model_parts:
                raise Exception(
                    "Duplicate sub_model_part_conditions entry for '{}'".format(sub_model_part_name))

            assigned_sub_model_parts.add(sub_model_part_name)

    def _UsesCustomConditionAssignments(self):
        return self.settings["sub_model_part_conditions"].size() > 0

    def _GetCustomConditionAssignments(self):
        custom_assignments = []
        assignments_settings = self.settings["sub_model_part_conditions"]

        for i in range(assignments_settings.size()):
            assignment = assignments_settings[i]
            custom_assignments.append((
                assignment["sub_model_part_name"].GetString(),
                assignment["condition_name"].GetString()))

        return custom_assignments

    def _PrepareThermalModelPartWithCustomConditions(self):
        thermal_solver = self.thermal_solver
        thermal_model_part = thermal_solver.main_model_part
        assign_neighbour_elements = thermal_solver.settings["assign_neighbour_elements_to_conditions"].GetBool()

        if not thermal_solver.is_restarted():
            materials_imported = thermal_solver.import_materials()
            if materials_imported:
                KratosMultiphysics.Logger.PrintInfo(
                    "::[ConvectionDiffusionSolver]:: ",
                    "Materials were successfully imported.")
            else:
                KratosMultiphysics.Logger.PrintInfo(
                    "::[ConvectionDiffusionSolver]:: ",
                    "Materials were not imported.")

            global_replace_settings = self._CreateThermalReplaceSettings(
                thermal_model_part,
                element_name=thermal_solver.settings["element_replace_settings"]["element_name"].GetString(),
                condition_name=thermal_solver.settings["element_replace_settings"]["condition_name"].GetString())
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(
                thermal_model_part,
                global_replace_settings).Execute()

            for sub_model_part_name, condition_name in self._GetCustomConditionAssignments():
                sub_model_part = self._GetThermalSubModelPart(sub_model_part_name)

                if sub_model_part.NumberOfConditions() == 0:
                    raise Exception(
                        "Configured thermal sub model part '{}' contains no conditions to replace.".format(
                            sub_model_part.FullName()))

                condition_replace_settings = self._CreateThermalReplaceSettings(
                    sub_model_part,
                    condition_name=condition_name)
                KratosMultiphysics.ReplaceElementsAndConditionsProcess(
                    sub_model_part,
                    condition_replace_settings).Execute()

                KratosMultiphysics.Logger.PrintInfo(
                    "::[CoupledFluidThermalSolverBoussinesq]::",
                    "Assigned '{}' to '{}'".format(
                        condition_replace_settings["condition_name"].GetString(),
                        sub_model_part.FullName()))

            tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
            throw_errors = False
            flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
            if assign_neighbour_elements:
                flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
            else:
                flags |= (tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS).AsFalse()
            tmoc(thermal_model_part, throw_errors, flags).Execute()

            thermal_solver._set_and_fill_buffer()

        if thermal_model_part.IsDistributed():
            thermal_solver.distributed_model_part_importer.CreateCommunicators()

        if thermal_solver.settings["echo_level"].GetInt() > 0:
            KratosMultiphysics.Logger.PrintInfo(self.model)

        KratosMultiphysics.Logger.PrintInfo(
            "::[ConvectionDiffusionSolver]::",
            "ModelPart prepared for Solver.")

    def _GetThermalSubModelPart(self, sub_model_part_name):
        thermal_model_part = self.thermal_solver.main_model_part
        full_sub_model_part_name = "{}.{}".format(thermal_model_part.FullName(), sub_model_part_name)

        if self.model.HasModelPart(full_sub_model_part_name):
            sub_model_part = self.model.GetModelPart(full_sub_model_part_name)
        elif thermal_model_part.HasSubModelPart(sub_model_part_name):
            sub_model_part = thermal_model_part.GetSubModelPart(sub_model_part_name)
        elif self.model.HasModelPart(sub_model_part_name):
            sub_model_part = self.model.GetModelPart(sub_model_part_name)
        else:
            raise Exception(
                "Thermal sub model part '{}' was not found in '{}'".format(
                    sub_model_part_name,
                    thermal_model_part.FullName()))

        if sub_model_part.GetRootModelPart().FullName() != thermal_model_part.GetRootModelPart().FullName():
            raise Exception(
                "Configured sub model part '{}' does not belong to thermal model part '{}'".format(
                    sub_model_part.FullName(),
                    thermal_model_part.FullName()))

        return sub_model_part

    def _CreateThermalReplaceSettings(self, model_part, element_name=None, condition_name=None):
        replace_settings = KratosMultiphysics.Parameters("{}")

        if element_name is not None:
            replace_settings.AddEmptyValue("element_name").SetString(
                self._ExpandThermalElementName(model_part, element_name))
        if condition_name is not None:
            replace_settings.AddEmptyValue("condition_name").SetString(
                self._ExpandThermalConditionName(model_part, condition_name))

        return replace_settings

    def _ExpandThermalElementName(self, model_part, element_name):
        if element_name == "":
            return element_name

        domain_size = self._GetThermalDomainSize()
        generic_element_names = [
            "EulerianConvDiff",
            "LaplacianElement",
            "MixedLaplacianElement",
            "AdjointHeatDiffusionElement",
            "QSConvectionDiffusionExplicit",
            "DConvectionDiffusionExplicit",
            "AxisymmetricEulerianConvectionDiffusion",
            "EulerianConvDiffShockCapturing"
        ]

        if element_name in generic_element_names:
            num_nodes = self._GetEntityNumberOfNodes(model_part.Elements, model_part, domain_size + 1)
            return f"{element_name}{domain_size}D{num_nodes}N"

        return element_name

    def _ExpandThermalConditionName(self, model_part, condition_name):
        if condition_name == "":
            return condition_name

        domain_size = self._GetThermalDomainSize()
        generic_condition_names = [
            "FluxCondition",
            "ThermalFace",
            "AxisymmetricThermalFace",
            "LineCondition",
            "ConsistentFluxBoundaryCondition"
        ]

        if condition_name in generic_condition_names:
            num_nodes = self._GetEntityNumberOfNodes(model_part.Conditions, model_part, domain_size)
            return f"{condition_name}{domain_size}D{num_nodes}N"

        return condition_name

    def _GetThermalDomainSize(self):
        domain_size = self.thermal_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2, 3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")
        return domain_size

    def _GetEntityNumberOfNodes(self, entity_container, model_part, default_num_nodes):
        num_nodes = 0
        if len(entity_container) > 0:
            for entity in entity_container:
                num_nodes = len(entity.GetNodes())
                break

        num_nodes = model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes)
        if not num_nodes:
            num_nodes = default_num_nodes

        return num_nodes

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

        # Create Quasi-Newton convergence accelerator if requested
        if self._coupling_strategy == "QuasiNewton":
            self._convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self._qn_settings)
            self._convergence_accelerator.Initialize()

        KratosMultiphysics.Logger.PrintInfo(
            "::[CoupledFluidThermalSolverBoussinesq]::",
            "Boussinesq process initialized.")

    def _NormalizeConvergenceFlag(self, is_converged):
        if is_converged is None:
            return True

        return bool(is_converged)

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
        # Store securely inside non-historical variable AUX_TEMPERATURE
        aux_phi_var = ConvectionDiffusionApplication.AUX_TEMPERATURE
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
            phi_var, aux_phi_var, self.thermal_solver.main_model_part, self.thermal_solver.main_model_part, 0)

        fluid_converged   = True
        thermal_converged = True
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
            fluid_converged = self._NormalizeConvergenceFlag(self.fluid_solver.SolveSolutionStep())

            # --- Step 3: solve Convection-Diffusion → φ^k ------------------
            # The C-D solver reads VELOCITY from the now-updated fluid model part.
            thermal_converged = self._NormalizeConvergenceFlag(self.thermal_solver.SolveSolutionStep())

            # --- Step 4/5: Relaxation/Quasi-Newton and residual --------
            # Compute relative residual BEFORE relaxation/QN so we see the raw prediction delta
            rel_res = ConvectionDiffusionApplication.BoussinesqCouplingUtilities.ComputeRelativeResidual(
                self.thermal_solver.main_model_part, phi_var, aux_phi_var)

            KratosMultiphysics.Logger.PrintInfo(
                "::[CoupledFluidThermalSolverBoussinesq]::",
                "Coupling it={0:3d} | |Δφ|/|φ| = {1:.3e}".format(it, rel_res))

            if rel_res <= self._coupling_tol:
                is_coupling_converged = True
                KratosMultiphysics.Logger.PrintInfo(
                    "::[CoupledFluidThermalSolverBoussinesq]::",
                    "Coupling converged in {0} iteration(s).".format(it))
                break

            if self._coupling_strategy == "QuasiNewton" and it > self._picard_iters_first:
                num_nodes = self.thermal_solver.main_model_part.NumberOfNodes()
                r_vec = KratosMultiphysics.Vector(num_nodes)
                x_vec = KratosMultiphysics.Vector(num_nodes)
                
                # Build Kratos Vectors for the accelerator using C++ utility
                ConvectionDiffusionApplication.BoussinesqCouplingUtilities.ComputeQuasiNewtonUpdateVectors(
                    self.thermal_solver.main_model_part, phi_var, aux_phi_var, x_vec, r_vec)
                
                # Apply Newton-Raphson update
                self._convergence_accelerator.UpdateSolution(r_vec, x_vec)
                
                # Write back converged/relaxed values using C++ utility
                ConvectionDiffusionApplication.BoussinesqCouplingUtilities.UpdateConvergenceVariables(
                    self.thermal_solver.main_model_part, phi_var, aux_phi_var, x_vec)
                
                self._convergence_accelerator.FinalizeNonLinearIteration()
                
            else:
                # Classic Picard with fixed relaxation
                if self._omega < 1.0:
                    ConvectionDiffusionApplication.BoussinesqCouplingUtilities.ApplyRelaxation(
                        self.thermal_solver.main_model_part, phi_var, aux_phi_var, self._omega)

                # φ^{k+1} becomes φ^k for the next iteration
                KratosMultiphysics.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                    phi_var, aux_phi_var, self.thermal_solver.main_model_part, self.thermal_solver.main_model_part, 0)

        self._boussinesq_process.ExecuteFinalizeSolutionStep()

        if not is_coupling_converged:
            KratosMultiphysics.Logger.PrintWarning(
                "::[CoupledFluidThermalSolverBoussinesq]::",
                "Coupling loop reached max iterations ({0}) without convergence. "
                "Last |Δφ|/|φ| = {1:.3e}".format(self._max_coupling_it, rel_res))

        return fluid_converged and thermal_converged and is_coupling_converged


