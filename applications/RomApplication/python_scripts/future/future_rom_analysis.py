import json
from pathlib import Path

import numpy

import KratosMultiphysics
from KratosMultiphysics.RomApplication.future.rom_decoders import FixedDofsDecoder, LinearDecoder
from KratosMultiphysics.RomApplication.future.rom_future_solver import FutureRomSolver


def CreateFutureRomAnalysisInstance(cls, global_model, parameters, future_rom_settings=None, rom_basis_output_folder="rom_data"):
    class FutureRomAnalysis(cls):
        """Makes a standard Kratos analysis stage solve its steps with the FutureRomSolver.

        The parent analysis stage sets up the model part, the DOFs and the processes. Its solver
        is not called in the solution steps.
        """

        def Initialize(self):
            super().Initialize()

            self.future_rom_solver = FutureRomSolver(self._GetSolver().GetComputingModelPart(), self._GetFutureRomSettings())
            self.future_rom_solver.Initialize()
            self.future_rom_solver.SetDecoder(self._CreateDecoder())
            if self.future_rom_solver.run_hrom:
                self._SetHRomEntities()
            self.rom_solutions = []

        def _SetHRomEntities(self):
            """HROM entities and weights from the files saved by the HROM training (one column of weights per set, e.g. per physics)."""
            folder = Path(rom_basis_output_folder)
            weights_set = self.future_rom_solver.hrom_weights_set
            def Weights(name, number_of_entities):
                weights = numpy.load(folder / f"HROM_{name}Weights.npy").reshape(number_of_entities, -1)
                return weights[:, weights_set]
            element_ids = numpy.load(folder / "HROM_ElementIds.npy")
            condition_ids = numpy.load(folder / "HROM_ConditionIds.npy")
            self.future_rom_solver.SetHRomEntities(element_ids, Weights("Element", len(element_ids)), condition_ids, Weights("Condition", len(condition_ids)))

        def _GetFutureRomSettings(self):
            """Settings of the FutureRomSolver. The nonlinear iterations not given take the values of the solver of the analysis
            (max_iteration and, with a displacement criterion, its tolerances), so that the ROM converges as the legacy one.
            A linear analysis ("analysis_type": "linear") is solved in one iteration, as the legacy linear strategy does.
            The scheme moves the mesh if the solver of the analysis does ("move_mesh_flag"), as required by the elements that
            use the current coordinates (e.g. total Lagrangian)."""
            settings = future_rom_settings.Clone() if future_rom_settings is not None else KratosMultiphysics.Parameters("{}")
            # The validated settings of the solver, as "analysis_type" and "move_mesh_flag" are usually taken from its defaults
            validated_solver_settings = self._GetSolver().settings
            if not settings.Has("max_iterations") and validated_solver_settings.Has("analysis_type") and validated_solver_settings["analysis_type"].GetString() == "linear":
                settings.AddInt("max_iterations", 1)

            solver_settings = self.project_parameters["solver_settings"]
            inherited = [("max_iterations", "max_iteration")]
            if solver_settings.Has("convergence_criterion") and solver_settings["convergence_criterion"].GetString() == "displacement_criterion":
                inherited += [("relative_tolerance", "displacement_relative_tolerance"), ("absolute_tolerance", "displacement_absolute_tolerance")]
            for name, legacy_name in inherited:
                if not settings.Has(name) and solver_settings.Has(legacy_name):
                    settings.AddValue(name, solver_settings[legacy_name])

            if validated_solver_settings.Has("move_mesh_flag"):
                if not settings.Has("scheme_settings"):
                    scheme_settings = FutureRomSolver.GetDefaultParameters()["scheme_settings"].Clone()
                    scheme_settings.RemoveValue("move_mesh")
                    settings.AddValue("scheme_settings", scheme_settings)
                if not settings["scheme_settings"].Has("move_mesh"):
                    settings["scheme_settings"].AddValue("move_mesh", validated_solver_settings["move_mesh_flag"])
            return settings

        def _CreateDecoder(self):
            """Linear decoder from the basis saved by the RomManager, taking the initial state as origin. With
            "fixed_dofs_decoding": "prescribed", wrapped in a FixedDofsDecoder."""
            folder = Path(rom_basis_output_folder)
            with open(folder / "RomParameters.json", 'r') as parameter_file:
                rom_settings = json.load(parameter_file)["rom_settings"]
            nodal_unknowns = rom_settings["nodal_unknowns"]
            node_ids = numpy.load(folder / "NodeIds.npy")
            right_modes = numpy.load(folder / "RightBasisMatrix.npy")
            if right_modes.ndim == 1:
                right_modes = right_modes.reshape(-1, 1)
            right_modes = right_modes[:, :rom_settings["number_of_rom_dofs"]]

            # The saved basis has one row per node and nodal unknown. Reorder them as the effective DOF set
            node_id_to_index = {node_id: index for index, node_id in enumerate(node_ids)}
            unknown_to_index = {name: index for index, name in enumerate(nodal_unknowns)}
            dof_set = self.future_rom_solver.GetEffectiveDofSet()
            rows = [node_id_to_index[dof.Id()] * len(nodal_unknowns) + unknown_to_index[dof.GetVariable().Name()] for dof in dof_set]

            decoder = LinearDecoder(right_modes[rows, :], numpy.array(dof_set.GetValues()))
            if self.future_rom_solver.fixed_dofs_decoding == "prescribed":
                # The fixed DOFs take their prescribed values instead of the ones of the basis
                decoder = FixedDofsDecoder(decoder, dof_set)
            return decoder

        def InitializeSolutionStep(self):
            self.PrintAnalysisStageProgressInformation()

            self.ApplyBoundaryConditions()
            self.ChangeMaterialProperties()
            self.future_rom_solver.InitializeSolutionStep()

        def SolveSolutionStep(self):
            is_converged = self.future_rom_solver.SolveSolutionStep()
            self.rom_solutions.append(self.future_rom_solver.GetRomSolution().copy())
            return is_converged

        def FinalizeSolutionStep(self):
            self.future_rom_solver.FinalizeSolutionStep()

            for process in self._GetListOfProcesses():
                process.ExecuteFinalizeSolutionStep()

    return FutureRomAnalysis(global_model, parameters)
