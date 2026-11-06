import json
from pathlib import Path

import numpy

import KratosMultiphysics
from KratosMultiphysics.RomApplication.future.rom_decoders import FixedDofsDecoder, LinearDecoder
from KratosMultiphysics.RomApplication.future.rom_local_decoder import LocalRomDecoder, NearestCentreSelector, ProcessInfoSelector
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
            self.rom_solutions = []

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
            "fixed_dofs_decoding": "prescribed", wrapped in a FixedDofsDecoder.

            With "local_roms" in the RomParameters.json, a LocalRomDecoder with one linear decoder per local model, each read from its
            own subfolder ("models") with the files of a single one (RightBasisMatrix.npy, and optionally ReferenceSnapshot.npy as its
            origin, else the initial state), and the "selector" of the model of each step:
            {"type": "process_info", "variable": "TIME", "thresholds": [...]} or {"type": "nearest_centre", "centres": "<file>.npy",
            "margin": 0.0}, with the centres saved as the bases (one row per node and nodal unknown, one column per model)."""
            folder = Path(rom_basis_output_folder)
            with open(folder / "RomParameters.json", 'r') as parameter_file:
                rom_parameters = json.load(parameter_file)
            rom_settings = rom_parameters["rom_settings"]
            rows = self._GetBasisRows(folder, rom_settings["nodal_unknowns"])
            u0 = numpy.array(self.future_rom_solver.GetEffectiveDofSet().GetValues())

            local_roms = rom_parameters.get("local_roms")
            if local_roms is None:
                return self._CreateLinearDecoder(folder, rom_settings["number_of_rom_dofs"], rows, u0, read_origin=False)

            models = local_roms["models"]
            number_of_rom_dofs = local_roms.get("number_of_rom_dofs", [rom_settings["number_of_rom_dofs"]] * len(models))
            decoders = [self._CreateLinearDecoder(folder / name, n, rows, u0, read_origin=True) for name, n in zip(models, number_of_rom_dofs)]
            return LocalRomDecoder(decoders, self._CreateModelSelector(local_roms["selector"], folder, rows))

        def _CreateLinearDecoder(self, folder, number_of_rom_dofs, rows, u0, read_origin):
            right_modes = numpy.load(folder / "RightBasisMatrix.npy")
            if right_modes.ndim == 1:
                right_modes = right_modes.reshape(-1, 1)
            right_modes = right_modes[:, :number_of_rom_dofs]
            origin = numpy.load(folder / "ReferenceSnapshot.npy")[rows] if read_origin and (folder / "ReferenceSnapshot.npy").exists() else u0
            decoder = LinearDecoder(right_modes[rows, :], origin)
            if self.future_rom_solver.fixed_dofs_decoding == "prescribed":
                # The fixed DOFs take their prescribed values instead of the ones of the basis
                decoder = FixedDofsDecoder(decoder, self.future_rom_solver.GetEffectiveDofSet())
            return decoder

        def _CreateModelSelector(self, settings, folder, rows):
            if settings["type"] == "process_info":
                return ProcessInfoSelector(self._GetSolver().GetComputingModelPart(), settings["variable"], settings["thresholds"])
            if settings["type"] == "nearest_centre":
                return NearestCentreSelector(numpy.load(folder / settings["centres"])[rows, :], settings.get("margin", 0.0))
            raise Exception("Unknown local ROM selector '{}'. Available options are 'process_info' and 'nearest_centre'.".format(settings["type"]))

        def _GetBasisRows(self, folder, nodal_unknowns):
            """The saved bases have one row per node and nodal unknown. Their rows in the order of the effective DOF set"""
            node_id_to_index = {node_id: index for index, node_id in enumerate(numpy.load(folder / "NodeIds.npy"))}
            unknown_to_index = {name: index for index, name in enumerate(nodal_unknowns)}
            dof_set = self.future_rom_solver.GetEffectiveDofSet()
            return [node_id_to_index[dof.Id()] * len(nodal_unknowns) + unknown_to_index[dof.GetVariable().Name()] for dof in dof_set]

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
