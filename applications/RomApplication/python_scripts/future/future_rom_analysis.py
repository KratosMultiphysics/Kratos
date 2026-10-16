import json
from pathlib import Path

import numpy

import KratosMultiphysics
from KratosMultiphysics.RomApplication.future.rom_decoders import LinearDecoder
from KratosMultiphysics.RomApplication.future.rom_future_solver import FutureRomSolver


def CreateFutureRomAnalysisInstance(cls, global_model, parameters, future_rom_settings=None, rom_basis_output_folder="rom_data"):
    class FutureRomAnalysis(cls):
        """Makes a standard Kratos analysis stage solve its steps with the FutureRomSolver.

        The parent analysis stage sets up the model part, the DOFs and the processes. Its solver
        is not called in the solution steps.
        """

        def Initialize(self):
            super().Initialize()

            self.future_rom_solver = FutureRomSolver(self._GetSolver().GetComputingModelPart(), future_rom_settings)
            self.future_rom_solver.Initialize()
            self.future_rom_solver.SetDecoder(self._CreateDecoder())
            self.rom_solutions = []

        def _CreateDecoder(self):
            """Linear decoder from the basis saved by the RomManager, taking the initial state as origin."""
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

            return LinearDecoder(right_modes[rows, :], numpy.array(dof_set.GetValues()))

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
