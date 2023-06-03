import importlib

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication import new_python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy


def CreateFomRomAnalysisInstance(cls, model, fom_project_parameters, rom_project_parameters):
    class FomRomAnalysis(cls):

        def __init__(self,model,fom_project_parameters,rom_project_parameters):
            self.rom_parameters = rom_project_parameters
            self.n_nodal_unknowns = len(self.rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray())
            super().__init__(model,fom_project_parameters)
            self._is_rom = False

        #### Internal functions ####
        def _CreateSolver(self):
            self._CreateRomSolver()
            self._CreateFomSolver()
            return self.fom_solver

        #### Internal functions ####
        def _CreateRomSolver(self):
            self.rom_solver_parameters = self.project_parameters.Clone()
            """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
            # Set the ROM settings in the "solver_settings" of the solver introducing the physics
            self.rom_solver_parameters["solver_settings"].AddValue("rom_settings", self.rom_parameters["rom_settings"])

            # ROM solving strategy
            self.solving_strategy = self.rom_parameters["projection_strategy"].GetString() if self.rom_parameters.Has("projection_strategy") else "galerkin"
            self.rom_solver_parameters["solver_settings"].AddString("projection_strategy",self.solving_strategy)

            # Add or remove parameters depending on the solving strategy
            ##LSPG
            if self.solving_strategy=="lspg":
                self.rom_solver_parameters["solver_settings"]["rom_settings"].AddBool("train_petrov_galerkin", self.train_petrov_galerkin)
            ##Petrov Galerkin
            if self.solving_strategy=="petrov_galerkin":
                self.petrov_galerkin_rom_dofs = self.rom_solver_parameters["solver_settings"]["rom_settings"]["petrov_galerkin_number_of_rom_dofs"].GetInt()
            else:
                self.rom_solver_parameters["solver_settings"]["rom_settings"].RemoveValue("petrov_galerkin_number_of_rom_dofs")

            # ROM assembling strategy
            self.assembling_strategy = self.rom_parameters["assembling_strategy"].GetString() if self.rom_parameters.Has("assembling_strategy") else "global"
            self.rom_solver_parameters["solver_settings"].AddString("assembling_strategy",self.assembling_strategy)

            # Create the ROM solver
            self.rom_solver = solver_wrapper.CreateSolver(self.model, self.rom_solver_parameters)

        def _CreateFomSolver(self):
            self.fom_solver = super()._CreateSolver()

        def _GetSimulationName(self):
            if self._is_rom:
                return "::[ROM Simulation]:: "
            else:
                return "::[FOM Simulation]" + super()._GetSimulationName()

        def _GetComputingModelPart(self):
            if self._is_rom:
                return self.rom_solver.GetComputingModelPart()
            else:
                return self.fom_solver.GetComputingModelPart()

        def SetSolverToRom(self):
            self._is_rom = True
            self._solver = self.rom_solver

        def SetSolverToFom(self):
            self._is_rom = False
            self._solver = self.fom_solver

        def UpdateRomBases(self,snapshots_matrix,svd_truncation_tolerance=1e-6):
            # Calculate the randomized SVD of the snapshots matrix
            u,_,_,eSVD = RandomizedSingularValueDecomposition().Calculate(snapshots_matrix, svd_truncation_tolerance)
            num_rom_dofs = numpy.shape(u)[1]

            node_counter = 0
            for node in self._GetComputingModelPart().Nodes:
                aux = KratosMultiphysics.Matrix(self.n_nodal_unknowns, num_rom_dofs)
                for j in range(self.n_nodal_unknowns):
                    for i in range(num_rom_dofs):
                        aux[j,i] = u.item((node_counter * self.n_nodal_unknowns +j,i))
                node.SetValue(KratosROM.ROM_BASIS, aux ) # ROM basis
                node_counter += 1

            self._GetComputingModelPart().ProcessInfo.SetValue(KratosROM.NUM_ROM_BASIS,num_rom_dofs)

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

        def Finalize(self):
            super().Finalize()
            if self._is_rom:
                self.fom_solver.Finalize()
            else:
                self.rom_solver.Finalize()

    return FomRomAnalysis(model, fom_project_parameters, rom_project_parameters)