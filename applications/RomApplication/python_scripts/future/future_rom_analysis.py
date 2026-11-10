from pathlib import Path

import numpy

import KratosMultiphysics
from KratosMultiphysics.RomApplication.hrom_training_utility import HRomTrainingUtility
from KratosMultiphysics.RomApplication.future.rom_decoders import CreateDecoder
from KratosMultiphysics.RomApplication.future.rom_future_solver import FutureRomSolver


def CreateFutureRomAnalysisInstance(cls, global_model, parameters, future_rom_settings=None, rom_basis_output_folder="rom_data"):
    class FutureRomAnalysis(cls):
        """Makes a standard Kratos analysis stage solve its steps with the FutureRomSolver.

        The parent analysis stage sets up the model part, the DOFs and the processes. Its solver
        is not called in the solution steps.

        With "train_hrom" the residual of each entity is projected onto the Jacobian of the decoder at the end of each step, and the
        HROM weights are computed from these snapshots when the simulation finishes (see _TrainHRom).

        With "save_rom_solutions" the reduced coordinates of each step are saved in the ROM folder (RomSolutions.npy). With
        "replay_rom_solutions" they are read from there and imposed in each step instead of solving it (FutureRomSolver.ImposeRomSolution):
        e.g. the HROM training of a solved run, or of another range of entities, without solving it again.
        """

        def Initialize(self):
            super().Initialize()

            self.future_rom_solver = FutureRomSolver(self._GetSolver().GetComputingModelPart(), self._GetFutureRomSettings())
            self.future_rom_solver.Initialize()
            self.future_rom_solver.SetDecoder(self._CreateDecoder())
            if self.future_rom_solver.run_hrom:
                self._SetHRomEntities()
            self.rom_solutions = []
            self.hrom_residual_snapshots = []
            if self.future_rom_solver.replay_rom_solutions:
                self.rom_solutions_to_replay = numpy.load(Path(rom_basis_output_folder) / "RomSolutions.npy")

        def _SetHRomEntities(self):
            """HROM entities and weights from the files saved by the HROM training (one column of weights per set, e.g. per physics)."""
            folder = Path(rom_basis_output_folder)
            weights_set = self.future_rom_solver.hrom_weights_set
            def Weights(name, number_of_entities):
                weights = numpy.load(folder / f"HROM_{name}Weights.npy")
                if weights.size == 0:
                    # No entity of this kind selected (reshape cannot infer the number of sets)
                    weights = weights.reshape(0, weights.shape[1] if weights.ndim > 1 else 1)
                else:
                    weights = weights.reshape(number_of_entities, -1)
                return weights[:, weights_set]
            element_ids = numpy.load(folder / "HROM_ElementIds.npy")
            condition_ids = numpy.load(folder / "HROM_ConditionIds.npy")
            self.future_rom_solver.SetHRomEntities(element_ids, Weights("Element", len(element_ids)), condition_ids, Weights("Condition", len(condition_ids)))

        def _TrainHRom(self):
            """HROM weights from the residual snapshots (ECM) and HROM model part, by the HRomTrainingUtility with the "hrom_settings"
            of the RomParameters.json. The weights are saved in the ROM folder as numpy files, the ones read by _SetHRomEntities."""
            folder = Path(rom_basis_output_folder)
            with open(folder / "RomParameters.json", 'r') as parameter_file:
                rom_parameters = KratosMultiphysics.Parameters(parameter_file.read())
            rom_parameters.AddString("rom_basis_output_name", "RomParameters")
            rom_parameters.AddString("rom_basis_output_folder", str(folder))
            if rom_parameters["hrom_settings"].Has("hrom_format") and rom_parameters["hrom_settings"]["hrom_format"].GetString() != "numpy":
                raise Exception("The HROM training of the Future ROM requires 'hrom_format' set to 'numpy'.")
            if self.project_parameters.Has("modelers"):
                rom_parameters["hrom_settings"].AddValue("modelers", self.project_parameters["modelers"])
            hrom_training_utility = HRomTrainingUtility(self._GetSolver(), rom_parameters)
            hrom_training_utility.time_step_residual_matrix_container = self.hrom_residual_snapshots
            hrom_training_utility.CalculateAndSaveHRomWeights()
            hrom_training_utility.CreateHRomModelParts()

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
            """The decoder saved in the ROM folder (see CreateDecoder in rom_decoders.py), with its rows in the order of the effective DOF set."""
            return CreateDecoder(rom_basis_output_folder, self.future_rom_solver.GetEffectiveDofSet())

        def InitializeSolutionStep(self):
            self.PrintAnalysisStageProgressInformation()

            self.ApplyBoundaryConditions()
            self.ChangeMaterialProperties()
            self.future_rom_solver.InitializeSolutionStep()

        def SolveSolutionStep(self):
            if self.future_rom_solver.replay_rom_solutions:
                step = len(self.rom_solutions)
                if step >= len(self.rom_solutions_to_replay):
                    raise Exception("The replayed RomSolutions.npy has {} steps, the simulation asks for step {}.".format(len(self.rom_solutions_to_replay), step + 1))
                self.future_rom_solver.ImposeRomSolution(self.rom_solutions_to_replay[step])
                self.rom_solutions.append(self.future_rom_solver.GetRomSolution().copy())
                return True
            is_converged = self.future_rom_solver.SolveSolutionStep()
            self.rom_solutions.append(self.future_rom_solver.GetRomSolution().copy())
            return is_converged

        def FinalizeSolutionStep(self):
            self.future_rom_solver.FinalizeSolutionStep()
            if self.future_rom_solver.train_hrom:
                self.hrom_residual_snapshots.append(self.future_rom_solver.GetProjectedEntityResiduals())

            for process in self._GetListOfProcesses():
                process.ExecuteFinalizeSolutionStep()

        def Finalize(self):
            super().Finalize()
            if self.future_rom_solver.save_rom_solutions:
                numpy.save(Path(rom_basis_output_folder) / "RomSolutions.npy", numpy.array(self.rom_solutions))
            if self.future_rom_solver.train_hrom:
                self._TrainHRom()

    return FutureRomAnalysis(global_model, parameters)
